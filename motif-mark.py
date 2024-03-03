#!/usr/bin/env python

import argparse
import cairo
import math
import re

# Global constants for image output
LEFT_MARGIN = 20
GENE_HEIGHT = 100
VERTICAL_BUFFER = GENE_HEIGHT / 5

###########
# Classes #
###########

class Gene:
    '''The gene sequence.'''
    def __init__(self):
        self.sequence = ""
        self.width = 0

    def add_sequence(self, seq):
        self.sequence += seq
        self.width = len(self.sequence)

    def draw_gene(self, context: cairo.Context, sequence_length):
        x = LEFT_MARGIN
        y = GENE_HEIGHT * sequence_length + (GENE_HEIGHT / 2)
        context.set_source_rgb(0,0,0)
        context.set_line_width(1)
        context.move_to(x,y)
        context.line_to(x + self.width, y)
        context.stroke()

class Drawing:
    '''A drawing sequence.'''
    def __init__ (self, gene):
        self.gene = gene
        self.context = context
    
    def draw(self, context):
        self.gene.draw_gene(context)
        
# class Exon:
#     '''This is how a exon is drawn.'''
#     def __init__(self, start: int, end:int, y_offset:int):

#         ## Data (Attributes) ##
#         self.start = start
#         self.end = end
#         self.y_offset = y_offset

#     def draw_exon(self, context: cairo.Context):
#         context = cairo.Context(surface)
#         context.set_source_rgb(0, 0, 0)
#         context.move_to(LEFT_MARGIN+self.start)
#         context.line_to(self.end)
#         context.stroke()
#         #surface.finish()
#         surface.write_to_png("plot.png")
#         my_exon = Exon(self, context: cairo.Context, beginning, end, exon_length)
#         my_exon.draw_exon(context)

########
# Main #
########
def get_args():
    parser= argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Input fasta", required=True)
    parser.add_argument("-m", "--motifs", help="Input motifs", required=True)
    return parser.parse_args()

args = get_args()
f=args.fasta
m=args.motifs 

def oneline_fasta(file:str,intermediate:str):
        '''This function takes a multi line fasta and converts it into a fasta that has a header line and one sequence line '''
        with open(file, "r")as fh, open(intermediate, "w") as fout:
            print(fh.readline(), end="", file=fout) #this makes it so the print statement wont add a newline character 
            for line in fh:
                line=line.strip('\n')
                if line.startswith('>'): #for every header ex for the first one 
                    print(f'\n{line}',file=fout)
                else:
                    print(line,end="",file=fout)
            print("",file=fout)

oneline_fasta(f, "intermediate.fasta")

# Initialize an empty dictionary
sequence_dict = {}

# Open the file "intermediate.fasta" and read its contents
with open("intermediate.fasta") as fasta_file:
    # Initialize variables to store the gene symbol and sequence
    gene_symbol = None
    sequence = ""
    
    # Iterate through each line in the file
    for line in fasta_file:
        # Strip any leading/trailing whitespace
        line = line.strip()
        
        # Check if the line starts with '>'
        if line.startswith('>'):
            # If there's a gene symbol and sequence, add them to the dictionary
            if gene_symbol and sequence:
                sequence_dict[gene_symbol] = sequence
                
            # Extract the new gene symbol
            gene_symbol = line.split()[0][1:]
            
            # Reset the sequence for the new gene
            sequence = ""
        else:
            # Concatenate the lines to form the sequence
            sequence += line
    
    # Add the last gene's sequence to the dictionary
    if gene_symbol and sequence:
        sequence_dict[gene_symbol] = sequence

# Print the swapped dictionary to verify
#print(sequence_dict)

for gene_symbol in sequence_dict:
    # Get the DNA sequence for the current gene symbol
    sequence = sequence_dict[gene_symbol]
    
    # Calculate the length of the sequence
    sequence_length = len(sequence)

    # Print the gene symbol and its corresponding sequence length
    #print(f"Gene: {gene_symbol}, Sequence Length: {sequence_length}")

    # Create a Cairo surface and context
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, 500, 500)
    context = cairo.Context(surface)

    # Set background color to white 
    context.set_source_rgb(1, 1, 1)  # White color
    context.paint()

    # Create a Gene object for the current sequence
    gene = Gene()
    gene.add_sequence(sequence)

    # Draw the gene
    gene.draw_gene(context, sequence_length)
    
surface.write_to_png("gene_sequence.png")


# def exon_finder(sequence):
#     '''Returns the start and stop positions of capitialized stretches within the sequence.'''
#     exon_sequence = re.finditer("[A-Z]+", sequence)
#     # returns start and stop position, 0 indexed
#     for start_end in exon_sequence:
#         exon_start_end = start_end.span()
#         return exon_start_end

# #make the exon object
#     exon = exon_finder(sequence)
#     exon_start = exon[0]
#     exon_end = exon[1]
#     exon_object =Exon(exon_start,exon_end,y_offset)



