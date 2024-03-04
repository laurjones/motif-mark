#!/usr/bin/env python

import argparse
import cairo
import math
import re

# Global constants for image output
LEFT_MARGIN = 20
GENE_HEIGHT = 100
VERTICAL_BUFFER = GENE_HEIGHT / 5
# Shifting genes downward 
Y_OFFSET = 50

###########
# Classes #
###########

class Gene:
    '''The gene sequence.'''
    def __init__(self, gene_number: int, gene_name: str):
        # self.sequence = ""
        self.width = 0
        self.gene_number = gene_number  # Initialize gene_number
        self.gene_name = gene_name
        # self.LEFT_MARGIN = 10  # Define LEFT_MARGIN
        # self.GENE_HEIGHT = 20 # Define GENE_HEIGHT
    
    def __str__(self) -> str:
        return f"My Gene: {self.width=}, {self.gene_number=}, {self.gene_name}"
    
    def add_sequence(self, seq):
        '''Function to determine the number of nucleotides in my sequence (width).'''
        # self.sequence = seq
        self.width = len(seq)

    def draw_gene(self, context: cairo.Context, gene_symbol):
        # self.gene_number += 1 # Increment gene_number
        x = LEFT_MARGIN 
        # y = GENE_HEIGHT * sequence_length + (GENE_HEIGHT / 2)
        y = GENE_HEIGHT * self.gene_number + Y_OFFSET
        context.set_source_rgb(0,0,0)
        context.set_line_width(1)
        context.move_to(x,y)
        context.show_text(gene_symbol)
        #context.show_text("jason rules")
        context.line_to(x + self.width, y)
        context.stroke()
        print(f'debug 123 {x=}, {y=}, {self.width=}')

class Exon:
    '''This is how a exon is drawn.'''
    def __init__(self, exon_start: int, exon_end:int, gene_number):
        ## Data (Attributes) ##
        self.exon_start = exon_start
        self.exon_end = exon_end
        self.gene_number = gene_number

    def draw_exon(self, context: cairo.Context):
        x = LEFT_MARGIN + self.exon_start
        y = GENE_HEIGHT * self.gene_number + Y_OFFSET
        context.set_source_rgb(0, 0.5, 0)
        context.set_line_width(4)
        context.move_to(x,y)
        context.line_to(LEFT_MARGIN + self.exon_end, y)
        context.stroke()
        #surface.finish()
        surface.write_to_png("plot.png")

class Drawing:
    '''A drawing sequence.'''
    def __init__ (self, gene, exon):
        self.gene = gene
        self.exon = exon
        # self.context = context
    
    def draw(self, context):
        self.gene.draw_gene(context)
        self.gene.draw_exon(context)

#############
# Functions #
#############
        
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

def exon_finder(sequence):
    '''Returns the start and stop positions of capitialized stretches within the sequence.'''
    exon_sequence = re.finditer("[A-Z]+", sequence)
    for start_end in exon_sequence:
        # eg exon_start_end = (5, 10)
        exon_start_end = start_end.span()
        return exon_start_end
#exon_finder("aaaaaTTTTTcccccGGGGGaaaaa")

#############
# Arguments #
#############
            
def get_args():
    parser= argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Input fasta", required=True)
    parser.add_argument("-m", "--motifs", help="Input motifs", required=True)
    return parser.parse_args()

args = get_args()
f=args.fasta
m=args.motifs 

########
# Main #
########

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

surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, 820, 500)
context = cairo.Context(surface)

# Set background color to white 
context.set_source_rgb(1, 1, 1)  # White color
context.paint()

# gene_number = 0
for gene_number, gene_symbol in enumerate(sequence_dict):
    print(f"debug 283 yes")
    # Get the DNA sequence for the current gene symbol
    sequence = sequence_dict[gene_symbol]
    
    # Calculate the length of the sequence
    sequence_length = len(sequence)

    # Print the gene symbol and its corresponding sequence length
    #print(f"Gene: {gene_symbol}, Sequence Length: {sequence_length}")

    # Create a Gene object for the current sequence
    gene = Gene(gene_number, gene_name = gene_symbol)
    gene.add_sequence(sequence)

    # Draw the gene
    print(f"debug 352 {str(gene)=}")
    gene.draw_gene(context, gene_symbol)
    # gene_number += 1

    exon = Exon(exon_start, exon_end, gene_number)
    #make the exon object
    exon_start_end = exon_finder(sequence)
    exon_start = exon_start_end[0]
    exon_end = exon_start_end[1]


surface.write_to_png("gene_sequence.png")




