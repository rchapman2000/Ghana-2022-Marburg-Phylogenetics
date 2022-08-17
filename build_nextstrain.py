import sys
import os
import re
import argparse as ap

from Bio import Entrez, SeqIO 

# Parses a date into the format YYYY-MM-DD and adds XXs
# to fields that are ambiguous.
#
# Input:
#   dt - an array containing date info
# 
# Output:
#   a string in the format YYYY-MM-DD
def parseDate(dt):
    dateArray = []
    
    # If the array has a length of 3, then year, month, and day were
    # provided.
    if len(dt) == 3:
        # No XX's need to be added
        dateArray = [dt[2],dt[1],dt[0]]
    # If the array has a length of 2, then only year andmonth were found.
    elif len(dt) == 2:
        ## Add XXs in the day position
        dateArray = [dt[1], dt[0], 'XX']
    # If the array only has a length of 1, then only the year was provided.
    else:
        # Add XXs in the day and month positions
        dateArray = [dt[0],'XX','XX']

    # There will not be a case where, the provided date array will 
    # be length 0, as the script skips records missing a date.

    # Join the array with dashes and return.
    return "-".join(dateArray)

def main():
    parser = ap.ArgumentParser()

    parser.add_argument('-i', '--input', required=True, type=str, \
        help='[Required] - input file name', \
        action='store', dest='infile') 
    parser.add_argument('-e', '--email', required = True, type=str, \
        help='[Required] - email address to be used to query NCBI', \
        action = 'store', dest='email')
    
    args = parser.parse_args()

    # Opens the input file
    i = open(args.infile, 'r')

    # Entrez requires an email, so the user's email
    # is provided.
    Entrez.email = args.email

    # Opens the ouput files: 
    # metadata.tsv - contains metadata for each sample including
    #                strainID, virus species, date, country, host
    #                host, accession number, authors, journal,
    #                and title
    #
    # sequences.fasta - contains the fasta sequences of each sample
    # countries.txt - a list of the unique countries that samples were
    #                 taken from
    metaout = open("metadata.tsv", "w+")
    seqout = open("sequences.fasta", "w+")
    countriesOut = open("countries.txt", "w+")

    # Defines an empty list to place countries of origin
    # in as samples are processed.
    countriesWritten = []

    # Writes the header of the metadata file.
    metaout.write("strain\tvirus\tdate\tcountry\thost\tAccession\tauthors\tjournal\ttitle\n")
    
    # Loops over each line in the input file
    for line in i:
        # Removes the newline character and splits the line
        # by the tab characters
        split = line.strip("\n").split("\t")
        print(split)

        # The header line starts with 'Strain Name', so
        # skip this line
        if split[0] != "Strain Name":
            
            # Checks to make sure that the date and country fields are
            # not missing. If they are, then the record is skipped.
            if split[5] != "-N/A-" and split[8] != '-N/A-':
                
                # Grabs fields from the Vipr metadata file
                virus = split[1]
                accession = split[2]
                # The date is supplied in the following formats
                # YYYY, MM/YYYY, or DD/MM/YYYY. This field is 
                # parsed to a list by splitting at the '/'
                # character and transformed into the YYYY-MM-DD format
                # including ambiguous charactera using the parseDate()
                # function
                date = parseDate(split[5].split("/"))
                host = split[6]
                country = split[8]
                
                # Defines variables for fields that will be
                # taken from the genbank record.
                authors = ''
                jounral = ''
                title = ''
                journal = ''

                # To ensure a unique strainID is applied to each sample,
                # the strainID is create from the accession number, stain 
                # supplied in the Vipr database, and date (separated by '|' characters)
                strain = "{0}|{1}|{2}".format(accession, split[0].replace(" ", ""), date)

                # Searches the nucleotide database for the given accesion number
                # and grabs the record from the handle
                recordHandle = Entrez.esearch(db="nucleotide", term="{0}[Accession]".format(accession))
                record = Entrez.read(recordHandle)

                # Fetches the genban record using the record ID found previously. Since an accession number
                # was searched, we should only get 1 results from the serach. Thus, we only need
                # to use the first record.
                idHandle = Entrez.efetch(db='nucleotide', id=record['IdList'][0], rettype='gb', retmode='txt')
                
                # Parses the genbank record using the SeqIO library
                seqRecord = SeqIO.read(idHandle, "genbank")
                print(accession)
                print(seqRecord)
                
                # Loops over the references in the records annotation
                for ref in seqRecord.annotations['references']:
                    # From the 'Director Submission' Reference,
                    # grabs the author, title, and journal.
                    if ref.title == 'Direct Submission':
                        print(ref)
                        author = ref.authors
                        title = ref.title
                        journal = ref.journal

                # Places the fields to be written in the output metadata file into
                # a list
                outfields = [strain, virus, date, country, host, accession, authors, journal, title]
                
                # Joins the outfields list with tabs, and writes the line to the output file
                metaout.write("\t".join(outfields) + "\n")

                # For the genome sequences, the header will need to be changed to the
                # strainID. To do this, we need to modify the id and description fields of
                # the genbank record. Thus, we make a copy of the record and change the fields
                # of the record.
                toWrite = seqRecord
                toWrite.id = strain
                toWrite.description = strain
                print(toWrite)
                
                # Writes the sequence in fasta format to sequences.fasta
                SeqIO.write(toWrite, seqout, 'fasta')

                # Checks to see if the sequence's country has been recorded
                if country not in countriesWritten:
                    # If not, the country is written to countries.txt
                    # and added to the list.
                    countriesOut.write(country + "\n")
                    countriesWritten.append(country)
    
    # Close file streams
    i.close()
    metaout.close()
    seqout.close()
    countriesOut.close()


if __name__ == "__main__":
    main()