#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
import os, sys, re

# Local imports
from utils import (
    Colors,
    err,
    fatal
)


def clean(s, remove=['"', "'"]):
    """Cleans a string to remove any defined leading or trailing characters.
    @param s <str>:
        String to clean.
    @param remove list[<str>]:
        List of characters to remove from beginning or end of string 's'.
    @return s <str>:
        Cleaned string
    """
    for c in remove:
        s = s.strip(c)
    return s


def check_file(file, ncols, delim='\t'):
    """Checks if a file is empty and has the expected number of columns.
    @param file <str>:
        Path to file to check file.
    @param ncols <int>:
        Number of expected columns, i.e. 2
    @param delim <str>:
        Delimiter of file, i.e. a tab
    @return header list[<str>]:
        Header of file represented as a list
    """
    # Check to see if the file is empty
    c = Colors()
    fh = open(file, 'r')
    try:
        header = [clean(col.lower().strip()) for col in next(fh).strip().split(delim)]
    except StopIteration:
        err('{}{}Error: --groups "{}" file is empty!{}'.format(c.bg_red, c.white, file, c.end))
        fatal('{}{}Please add Sample and Group information to the file and try again.{}'.format(c.bg_red, c.white, c.end))
    finally:
        fh.close()

    # Check for expected number of columns
    # on each line, report all problematic
    # lines and then exit 1 if any errors
    errors = False
    with open(file, 'r') as fh:
        linenumber = 0
        for line in fh:
            linenumber += 1
            linelist = [l.strip() for l in line.split(delim) if l.strip()]
            if len(linelist) == ncols or len(linelist) == 0:
                # Skip over lines with the expected
                # number of columns and empty line
                continue
            # Line is missing a column or it is not
            # tab delimited
            errors = True
            err(
                '{}{}Error: --groups "{}" file does not contain the expected number of columns at line {}! {}'.format(
                   c.bg_red, c.white, file, linenumber, c.end
                )
            )
            err('{}{}  └── Bad line contents: "{}" {}'.format(c.bg_red, c.white, line.rstrip(), c.end))
    if errors:
        fatal(
            '{}{}Fatal: Please correct these lines and check your file is tab-delimited! {}'.format(
                c.bg_red, c.white, c.end
            )
        )

    return header


def index(file, delim='\t', required = ['Sample', 'Group']):
    """Returns the index of expected columns in provided file. The groups 
    file is expected to have the following required columns. 
    @Required columns:
        - Sample, Group
    @param file <str>:
        Path to groups TSV file.
    @return tuple(indices <dict[int/None]>, hasHeader <boolean>):
        [0] Dictionary containing information the index of each required/optional column
        [1] Boolean to indicate whether file has a header
    """
    c = Colors()
    indices = {}
    has_header = True  
    
    # Check to see if the file is empty,
    # and it has the expected number of 
    # columns, i.e 2 
    header = check_file(file = file, ncols = len(required), delim = delim)
    # Parse the header to get the index of required fields
    # Get index of Sample, Group
    # columns for parsing the file
    warning = False 
    for col in required:
        try: indices[col] = header.index(col)
        except ValueError:
            warning = True 
            # Missing column names or header in groups file
            # This can also occur if the file is not actually 
            # a tab delimited file.
            # TODO: Add a check to see if the file is actually
            # a tab delimited file, i.e. a TSV file.
            has_header = False
            err('{}{}Warning: {} is missing the following column name: {} {}'.format(
                    c.bg_yellow, c.black, file, col, c.end
                )
            )
    if warning:
        err('{}{}  └── Making assumptions about columns in the groups file... 1=Sample, 2=Group {}'.format(
                c.bg_yellow, c.black,c.end
            )
        )
        # Setting column indexes to the following defaults:
        # 0 = Sample basename column
        # 1 = Group information column
        for i in range(len(required)):
            indices[required[i]] = i

    return indices, has_header


def groups(file, delim='\t'):
    """Reads and parses a sample sheet, groups.tsv, into a dictionary. 
    This file acts as a sample sheet to gather sample metadata and define 
    relationship betweeen groups of samples. This tab delimited file 
    contains two columns. One column for the basename of the sample and 
    lastly, one column for the name of the sample's group. It is worth
    noting that a sample can belong to more than one group. A 1:M sample
    to group relationship can be denoted by seperating muliptle groups
    with commas (i.e. ','). This group information is used downstream
    in the pipeline for differential methylation rules. Comparisons
    between groups can be made with a constrast.tsv file. This function
    returns a tuple containing a dictionary containing group to sample list.
    @param file <str>:
        Path to groups TSV file.
    @return pairs <dict[str]>:
        Dictionary containing ChIP-input pairs, where each key is ChIP 
        sample and its value is its matched input sample
    @return groups <dict[str]>:
        Dictionary containing group to samples, where each key is group 
        and its value is a list of samples belonging to that group
    @Example: group.tsv
        Sample        Group
        Sample_1	G1,G4
        Sample_2	G1,G4
        Sample_3	G1,G4
        Sample_4	G2,G5
        Sample_5	G2,G5
        Sample_6	G2,G5
        Sample_7	G3,G5
        Sample_8	G3,G5
        Sample_9	G3
        Sample_0	G3
    # Example data structure that is return given the file above
    >> groups
    {
        'G1': ['Sample_1', 'Sample_2', 'Sample_3'],
        'G2': ['Sample_4', 'Sample_5', 'Sample_6'],
        'G3': ['Sample_7', 'Sample_8', 'Sample_9', 'Sample_0'],
        'G4': ['Sample_1', 'Sample_2', 'Sample_3'],
        'G5': ['Sample_4', 'Sample_5', 'Sample_6', 'Sample_7', 'Sample_8'],
    }
    """
    # Get index of each required and 
    # optional column and determine if
    # the file has a header 
    indices, header = index(file)
    s_index = indices['Sample']
    g_index = indices['Group']
    # Parse sample and group 
    # information 
    groups = {}
    with open(file, 'r') as fh:
        if header:
            # Skip over header
            tmp = next(fh)
        for line in fh:
            linelist = [l.strip() for l in line.split(delim)]
            # Parse Sample information
            try:
                sample_name = clean(linelist[s_index])
                if not sample_name: continue  # skip over empty string
            except IndexError:
                continue  # No sample information, skip over line
            # Parse Group Information
            try:
                group = linelist[g_index]
                if not group: continue # skip over empty string
            except IndexError:
                continue  # No group information, skip over line
            # Check for multiple groups,
            # split on comma or semicolon
            multiple_groups = re.split(';|,',group)
            multiple_groups = [clean(g.strip()) for g in multiple_groups]
            for g in multiple_groups:
                if g not in groups:
                    groups[g] = []
                if sample_name not in groups[g]:
                    groups[g].append(sample_name)
    return groups


def contrasts(file, groups, delim='\t'):
    """Reads and parses the group comparison file, contrasts.tsv, into a 
    dictionary. This file acts as a config file to setup contrasts between
    two groups, where groups of samples are defined in the groups.tsv file.
    This information is used in differential analysis, like differential binding
    analysis or differential gene expression, etc. 
    @Example: contrasts.tsv
        G2  G1
        G4  G3
        G5  G1
    >> contrasts = contrasts('contrasts.tsv', groups = ['G1', 'G2', 'G3', 'G4', 'G5'])
    >> contrasts
    [
        ["G2",  "G1"],
        ["G4",  "G3"],
        ["G5",  "G1"]
    ]
    @param file <str>:
        Path to contrasts TSV file.
    @param groups list[<str>]:
        List of groups defined in the groups file, enforces groups exist.
    @return comparisons <list[list[str, str]]>:
        Nested list contain comparsions of interest.  
    """
    c = Colors()
    errors = []
    comparsions = []
    line_number = 0
    with open(file) as fh:
        for line in fh:
            line_number += 1
            linelist = [clean(l.strip()) for l in line.split(delim)]
            try:
                g1 = linelist[0]
                g2 = linelist[1]
                if not g1 or not g2: continue # skip over empty lines
            except IndexError:
                # Missing a group, need two groups to tango
                # This can happen if the file is NOT a TSV file,
                # and it is seperated by white spaces, :(  
                err(
                '{}{}Warning: {} is missing at least one group on line {}: {}{}'.format(
                    c.bg_yellow,
                    c.black,
                    file,
                    line_number,
                    line.strip(),
                    c.end
                    )
                )
                err('{}{}  └── Skipping over line, check if line is tab seperated... {}'.format(
                    c.bg_yellow,
                    c.black,
                    c.end)
                )
                continue
            # Check to see if groups where defined already,
            # avoids user errors and spelling errors
            for g in [g1, g2]:
                if g not in groups:
                    # Collect all error and report them at end
                    errors.append(g)
            
            # Add comparsion to list of comparisons
            if [g1, g2] not in comparsions:
                comparsions.append([g1, g2])

    if errors:    
        # One of the groups is not defined in groups
        err('{}{}Error: the following group(s) in "{}" are not defined in groups file! {}'.format(
            c.bg_red, 
            c.white,
            file,
            c.end)
        )
        fatal('{}{}  └── {} {}'.format(
            c.bg_red,
            c.white,
            ','.join(errors),
            c.end)
        )
    
    return comparsions


if __name__ == '__main__':
    # Testing groups TSV parser
    print('Parsing groups file...')
    grp2sample = groups(sys.argv[1])
    print(grp2sample)
    print('Parsing contrasts file...')
    comparsions = contrasts(sys.argv[2], groups=grp2sample.keys())
    print(comparsions)
