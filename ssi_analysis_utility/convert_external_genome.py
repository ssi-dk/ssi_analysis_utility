#!/usr/bin/env python3

__author__ = "David Smith"
__version__ = "1.0.0"
__email__ = "dsmith@tgen.org"
__modified_for_SSI_use_by = "Thor Bech Johannesen, thej@ssi.dk"

import logging
class GenomeStatus(object):
    """
    Contains and manipulates any generic data that is per-contig-position.
    This could be any single type of data, like actual bases, filter data,
    depth information, insertion data, pileups, etc.
    In the perl version, this object was originally a hash of lists, and
    converted to a hash of strings for performance.  In this version, it was
    originally a dictionary of strings, and converted to a dictionary of
    lists for performance and flexibility.
    Whether single characters, numbers, boolean values, or a mixture of data
    types are used does not seem to affect memory and performance.
    Storing lists per-contig-position with this class is a complicated
    affair, as several of the manipulation functions assume you mean to
    manipulate a continuous range of positions instead of a single position
    when you do that.
    """

    """
    The conventions used for what data is stored are as follows:
    Genomes:
        A, C, G, T, U, a, c, g, t, u:  The respective call.
        N, n:  Called "N" according to upstream analysis tools.
        X:  Not called by upstream analysis tools.
        . or empty string:  A deletion relative to reference.
        String of length >1:  An insertion relative to reference.
        Any other single letter:  A degeneracy.
        !:  An algorithmic failure.
    Duplicate region data:
        0:  Position not in a region that is duplicated within the reference.
        1:  Position is in a region that is duplicated.
        -:  Duplicate checking at this position was skipped by the user.
        !:  An algorithmic failure.
    Filters:
        Y:  This position passed its filter.
        N:  This position failed its filter.
        ?:  The filter could not be checked, and so the position is assumed
            to have failed.
        -:  The filter was not applicable, or skipped, or could not be checked
            for a known reason, and so is assumed to have passed.
        !:  An algorithmic failure.
    """

    # Arrays are zero-indexed, genome positions are one-indexed. Off-by-one errors? Never heard of 'em.
    def __init__(self):
        """
        Attributes:
            _status_data: The dictionary of lists that stores the actual genome.
            data. The keys of the dictionary are the contig names.  The lists
            correspond to the position data on that contig.  Genome position is
            list position + 1.
            _current_contig: Tracks the most recently-referenced contig, for
            convenience EG reading in fastas line-by-line.
        """
        self._status_data = {}
        self._current_contig = None

    # NOTE(jtravis): Does not raise exception if contig name is None or the empty string as documented
    def add_contig(self, contig_name):
        """
        Defines a new empty contig in the genome.
        By default, if an unrecognized contig is encountered, a new empty
        contig will be created and then acted upon.
        Otherwise, add_contig must be called on a new contig first, or an
        InvalidContigName will be thrown.

        Args:
            contig_name (str): Unique contig description.

        Raises:
            InvalidContigName: If contig_name is undefined.
        """
        if contig_name not in self._status_data:
            self._status_data[contig_name] = []
        self._current_contig = contig_name

    # NOTE(jtravis): unused parameter create_contig
    def set_current_contig(self, contig_name, create_contig=True):
        """
        Sets the most-recently-referenced contig without actually performing
        any action on the data.
        Can be called to return the current contig without changing it if
        given a contig_name of None.
        Will create the contig if it has not been encountered yet by
        default, or throw an InvalidContigName otherwise.

        Args:
            contig_name (str): Unique contig description or None to query the current contig name.
            create_contig (bool): If True and the contig does not exist, an empty contig will be created.

        Returns:
            str: Name of the last accessed contig or None.

        Raises:
            InvalidContigName: If create_contig is False and the contig does not exist.
        """
        if contig_name is None:
            contig_name = self._current_contig
        elif contig_name in self._status_data:
            self._current_contig = contig_name
        elif create_contig:
            self.add_contig(contig_name)
        else:
            raise InvalidContigName(contig_name, self.get_contigs())
        return contig_name

    def get_contigs(self):
        """
        Returns:
            list: Sorted list of contig names.
        """
        return sorted(self._status_data.keys())

    def append_contig(self, genome_data, contig_name=None):
        """
        Places the passed-in data at the position following the last
        defined position on the contig.  If passed a list, will give each
        item in the list its own position.

        Args:
            genome_data (list): List of nucleotide symbols.
            contig_name (str): Unique contig description.
        """
        contig_name = self.set_current_contig(contig_name)
        self._status_data[contig_name].extend(genome_data)

    def extend_contig(self, new_length, missing_range_filler, contig_name=None):
        """
        Ensures the contig is at least new_length positions long

        Args:
            new_length (int): Minimum contig length.
            missing_range_filler (str): Placeholder character for undefined areas at the end of the contig.
            contig_name (str): Unique contig description.
        """
        contig_name = self.set_current_contig(contig_name)
        if len(self._status_data[contig_name]) < new_length:
            self._status_data[contig_name].extend(
                [missing_range_filler] * ( new_length - len(self._status_data[contig_name]) ))

    def set_value(self, new_data, position_number, missing_range_filler="!", contig_name=None):
        """
        Sets the value at position_number on the contig.
        If passed a list, will change the continuous range of positions
        starting at position_number, one position per list item.
        Will extend the contig with missing_range_filler filling undefined
        values if the position to set is beyond the end of the contig.

        Args:
            new_data (str or list): Single or list of nucleotide symbols.
            position_number (int): 1-indexed contig position number.
            missing_range_filler (str): Filler for undefined regions before the set value. Modifies the data.
            contig_name (str): Unique contig description
        """
        contig_name = self.set_current_contig(contig_name)
        self.extend_contig(position_number, missing_range_filler, contig_name)
        if len(new_data) > 1:
            self._status_data[contig_name][position_number - 1:position_number - 1 + len(new_data)] = new_data
        else:
            self._status_data[contig_name][position_number - 1] = new_data

    def get_value(self, first_position, last_position=None, contig_name=None, filler_value=None):
        """
        Args:
            contig_name (str): Unique contig description.
            first_position (int): 1-indexed first position number.
            last_position (int): Optional last position to select a range or -1 to specify the end of the contig.
            filler_value (str): Optional filler for undefined regions beyond the genome data. Does not modify the data.

        Returns:
            Returns the nucleotide at first_position, list of values from
            first_position to last_position inclusive, or None.
        """
        contig_name = self.set_current_contig(contig_name)
        queried_value = filler_value
        if last_position is None:
            if first_position <= len(self._status_data[contig_name]):
                queried_value = self._status_data[contig_name][first_position - 1]
        else:
            queried_value = []
            if last_position == -1:
                last_position = len(self._status_data[contig_name])
            if last_position >= first_position and first_position <= len(self._status_data[contig_name]):
                queried_value = self._status_data[contig_name][first_position - 1:last_position]
                if filler_value is not None and len(queried_value) < last_position - first_position + 1:
                    queried_value.extend([filler_value] * ( last_position - first_position + 1 - len(queried_value) ))
        return queried_value

    def get_contig_length(self, contig_name=None):
        """
        Args:
            contig_name (str): Unique contig description.

        Returns:
            int: Number of positions defined in the contig
        """
        contig_name = self.set_current_contig(contig_name)
        return len(self._status_data[contig_name])

    def send_to_fasta_handle(self, output_handle, contig_prefix="", max_chars_per_line=80):
        """
        Assumes the genome data is in string format or stringifiable and one
        character per position, and then writes it in to the handle open for
        writing.  The file format is like a typical fasta were the genome
        data to be base calls (but no checks are performed).

        The contigs are sorted by name, not the order they were created.

        Args:
            output_handle (file object): File to append FASTA string.
            contig_prefix (str): Prefix for all contig names.
            max_chars_per_line (int): A positive value will limit the max chars per line.
        """
        for current_contig in self.get_contigs():
            output_handle.write(">" + contig_prefix + current_contig + "\n")
            if max_chars_per_line > 0:
                i = 0
                while ( max_chars_per_line * i ) < len(self._status_data[current_contig]):
                    output_handle.write(''.join(self._status_data[current_contig][
                                                ( max_chars_per_line * i ):( max_chars_per_line * ( i + 1 ) )]) + "\n")
                    i += 1
            else:
                output_handle.write(''.join(self._status_data[current_contig]) + "\n")

    def write_to_fasta_file(self, output_filename, contig_prefix="", max_chars_per_line=80):
        """
        Opens the passed filename and passes to send_to_fasta_handle.
        This is a separate function so that unit testing is easier, and
        file names or open file handles can be used as destinations.

        Args:
            output_filename (str): Output filename.
            contig_prefix (str): Prefix for all contig names.
            max_chars_per_line (int): A positive value will limit the max contig chars per line.
        """
        with open(output_filename, 'w') as output_handle:
            self.send_to_fasta_handle(output_handle, contig_prefix, max_chars_per_line)

class Genome(GenomeStatus):
    """
    A special type of GenomeStatus where the genome information being stored
    is always actual base calls, as strings.
    """

    def __init__(self):
        """
        Attributes:
            _genome is an alias of _status_data, provided for code clarity
            when working with an actual genome
        """
        GenomeStatus.__init__(self)
        self._genome = self._status_data

    def set_call(self, new_data, first_position, missing_range_filler="X", contig_name=None):
        """ Alias of set_value, for code clarity """
        self.set_value(new_data, first_position, missing_range_filler, contig_name)

    def get_call(self, first_position, last_position=None, contig_name=None, filler_value="X"):
        """ Alias of get_value, for code clarity """
        return self.get_value( first_position, last_position, contig_name, filler_value)

    # NOTE: contig_prefix is unused
    def _import_fasta_line(self, line_from_fasta, contig_prefix=""):
        """
        Assumes the string passed in is a line from a fasta file, and
        populates the genome with the information contained.  Not meant
        to be called on any data except a full fasta file in order by
        line top to bottom.

        Args:
            line_from_fasta (str): the current line to parse
            contig_prefix (str): the prefix will be removed from the parsed contig name
        """
        import re

        # Parse the contig name discarding the prefix and surrounding whitespace characters
        contig_match = re.match(r'^>' + re.escape(contig_prefix) + r'([^\s]+)(?:\s|$)', line_from_fasta)
        if contig_match:
            self.add_contig(contig_match.group(1))
        else:
            # Parse the contig sequence discarding trailing whitespace characters
            data_match = re.match(r'^([A-Za-z.-]+)\s*$', line_from_fasta)
            if data_match:
                self.append_contig(list(data_match.group(1)))

    # contig_prefix is used by vcf_to_matrix to discard the frankenfasta contig name prefix.
    def import_fasta_file(self, fasta_filename, contig_prefix=""):
        """ Read in a fasta file.

        Args:
            fasta_filename (str): fasta file to import
            contig_prefix (str): the prefix will be removed from the parsed contig names
        """
        with open(fasta_filename, 'r') as fasta_handle:
            for line_from_fasta in fasta_handle:
                self._import_fasta_line(line_from_fasta, contig_prefix)

    @staticmethod
    def reverse_complement(dna_string):
        """
        Args:
            dna_string (str): nucleotide sequence to reverse complement

        Returns:
            string: nucleotide sequence reverse complement
        """
        return dna_string.translate(
            ''.maketrans('ABCDGHMNRSTUVWXYabcdghmnrstuvwxy', 'TVGHCDKNYSAABWXRtvghcdknysaabwxr'))[::-1]

    @staticmethod
    def simple_call(dna_string, allow_x=False, allow_del=False):
        """
        Standardizes the DNA call assumed to be the base at position one.
        Discards insertion data, changes 'U' to 'T', and changes degeneracies to 'N'.
        'X' and deletes are changed to 'N' by default.

        Args:
            dna_string (str): only the first position is considered
            allow_x (bool):
            allow_del (bool):

        Returns:
            string: 'A', 'C', 'G', 'T', or 'N' with optional 'X' and '.'
        """
        simple_base = 'N'
        if len(dna_string) > 0:
            simple_base = dna_string[0].upper()
        elif allow_del:
            simple_base = '.'
        if simple_base == 'U':
            simple_base = 'T'
        if simple_base not in ['A', 'C', 'G', 'T', 'X', '.']:
            simple_base = 'N'
        if not allow_x and ( simple_base == 'X' ):
            simple_base = 'N'
        if not allow_del and ( simple_base == '.' ):
            simple_base = 'N'
        return simple_base


def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Meant to be called from the pipeline automatically.")
    parser.add_argument("--nucmerpath", default="nucmer", help="Path to the 'nucmer' executable.")
    parser.add_argument("--nucmerargs", default="", help="Optional arguments to pass to the 'nucmer' executable.")
    parser.add_argument("--deltafilterpath", default="delta-filter", help="Path to the 'delta-filter' executable.")
    parser.add_argument("--deltafilterargs", default="", help="Optional arguments to pass to the 'delta-filter' executable.")
    parser.add_argument("--reference", required=True, help="Path to the reference fasta file.")
    parser.add_argument("--external", required=True, help="Path to the external genome fasta file.")
    parser.add_argument("--name", default="", help="Name of this external genome.")
    parser.add_argument("--outputdir", default="", help="Path to outputfolder containing delta and frankenfasta files")
    return parser.parse_args()


# This should eventually be moved to the main job manager section


def _update_genome_from_delta_data(franken_genome, external_genome, parser_state, distance_covered, is_external_insert):

    if distance_covered == -1:
        distance_covered = parser_state['final_pos'] - parser_state['reference_pos'] + 1
        is_external_insert = True
    if distance_covered > 0:
        if parser_state['external_is_reversed']:
            matching_segment = Genome.reverse_complement(''.join(
                external_genome.get_call(( parser_state['external_pos'] - distance_covered + 1 ),
                                         parser_state['external_pos'])))
        else:
            matching_segment = ''.join(external_genome.get_call(parser_state['external_pos'], (
                parser_state['external_pos'] + distance_covered - 1 )))
        franken_genome.set_call(matching_segment, parser_state['reference_pos'], 'X')
    parser_state['reference_pos'] = parser_state['reference_pos'] + distance_covered
    parser_state['external_pos'] = parser_state['external_pos'] + (
        -distance_covered if parser_state['external_is_reversed'] else distance_covered )
    if is_external_insert:
        parser_state['external_pos'] += -1 if parser_state['external_is_reversed'] else 1
    else:
        franken_genome.set_call('.', parser_state['reference_pos'], '!')
        parser_state['reference_pos'] += 1
    return parser_state


def _parse_delta_line(line_from_delta_file, franken_genome, external_genome, parser_state):
    import re

    line_match = re.match(r'^>([^ ]+) ([^ ]+) (\d+) \d+\s*$', line_from_delta_file)
    if line_match:
        current_contig = line_match.group(1)
        external_genome.set_current_contig(line_match.group(2))
        parser_state['contig_sizes'][current_contig] = int(line_match.group(3))
        franken_genome.add_contig(current_contig)
    else:
        line_match = re.match(r'^(\d+) (\d+) (\d+) (\d+) \d+ \d+ \d+\s*$', line_from_delta_file)
        if line_match:
            parser_state['reference_pos'] = int(line_match.group(1))
            parser_state['final_pos'] = int(line_match.group(2))
            parser_state['external_pos'] = int(line_match.group(3))
            parser_state['external_is_reversed'] = (
                True if ( parser_state['external_pos'] > int(line_match.group(4)) ) else False )
        else:
            line_match = re.match(r'^(\-?)(\d+)\s*$', line_from_delta_file)
            if line_match:
                distance_covered = int(line_match.group(2)) - 1
                is_external_insert = ( True if ( line_match.group(1) == '-' ) else False )
                parser_state = _update_genome_from_delta_data(franken_genome, external_genome, parser_state,
                                                              distance_covered, is_external_insert)
    return parser_state


def parse_delta_file(delta_filename, franken_genome, external_genome):
    parser_state = dict(zip(['contig_sizes', 'reference_pos', 'external_pos', 'final_pos', 'external_is_reversed'],
                            [dict(), None, None, None, None]))
    delta_handle = open(delta_filename, 'r')
    for line_from_delta_file in delta_handle:
        parser_state = _parse_delta_line(line_from_delta_file, franken_genome, external_genome, parser_state)
    delta_handle.close()
    for current_contig in franken_genome.get_contigs():
        franken_genome.extend_contig(parser_state['contig_sizes'][current_contig], 'X', current_contig)


def main():
    import subprocess
    import os

    commandline_args = _parse_args()

    external_genome = Genome()
    external_nickname = commandline_args.name
    external_genome.import_fasta_file(commandline_args.external)
    generate_delta_file(commandline_args.nucmerpath, commandline_args.nucmerargs, commandline_args.deltafilterpath,
                        commandline_args.deltafilterargs, external_nickname, commandline_args.reference,
                        commandline_args.external, commandline_args.outputdir)
    franken_genome = Genome()
    parse_delta_file(( os.path.join(commandline_args.outputdir,external_nickname + ".filtered.delta" )), franken_genome, external_genome)
    franken_genome.write_to_fasta_file(( os.path.join(commandline_args.outputdir,external_nickname + ".frankenfasta")), external_nickname + " ref:")

    #import pkg_resources
    #nasptool_path = pkg_resources.resource_filename('nasp', 'nasptool_linux_64')
    #with open(external_nickname + ".frankenfasta", 'w') as handle:
    #    return_code = subprocess.call([nasptool_path, "frankenfasta", external_nickname + ".filtered.delta"], stdout=handle)


if __name__ == "__main__":
    main()

