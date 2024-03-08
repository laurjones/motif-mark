#!/usr/bin/env python

import argparse
import cairo
import re

# Global constants for image output
LEFT_MARGIN = 20
GENE_HEIGHT = 100
VERTICAL_BUFFER = GENE_HEIGHT / 5
# Shifting genes downward 
Y_OFFSET = 50

# assign one color per motif
# key = motif (ygcy), value = color (0.9, 0.1, 0.1)
MOTIF_COLOR_DICT: dict[str, tuple[float, float, float]] = {} 
COLOR_PALETTE = [
        (0.9, 0.1, 0.1),  # Red
        (0.1, 0.9, 0.1),  # Green
        (0.1, 0.1, 0.9),  # Blue
        (0.9, 0.9, 0.1),   # Yellow
        (0.8, 0.3, 0.4),   # Mystery
        # Add more colors as needed
]

###########
# Classes #
###########

class Gene:
    '''This is how a gene sequence is drawn.'''
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
        context.move_to(x,y-35)
        context.show_text(gene_symbol)
        #context.show_text("jason rules")
        context.move_to(x,y)
        context.line_to(x + self.width, y)
        context.stroke()
        #print(f'debug 123 {x=}, {y=}, {self.width=}')

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
        context.set_source_rgb(0, 0, 0)
        context.set_line_width(15)
        context.move_to(x,y)
        context.line_to(LEFT_MARGIN + self.exon_end, y)
        context.stroke()
        #surface.finish()

class Motif:
    '''This is how a motif is drawn.'''
    def __init__(self, motif_start: int, motif_end:int, gene_number: int, motif_str: str):
        ## Data (Attributes) ##
        self.motif_start = motif_start
        self.motif_end = motif_end
        self.gene_number = gene_number
        self.motif_str = motif_str      # eg 'ygcy'
    
        self.motif_color = MOTIF_COLOR_DICT[motif_str]

    def __repr__(self) -> str:
        return f"Motif({self.motif_start}, {self.motif_end}, {self.gene_number}, '{self.motif_str}')"

    def draw_motif(self, context: cairo.Context):
        x = LEFT_MARGIN + self.motif_start
        y = GENE_HEIGHT * self.gene_number + Y_OFFSET
        context.set_source_rgb(*self.motif_color)  # Assuming motif_color is a tuple of RGB values
        context.set_line_width(8)
        context.move_to(x, y)
        context.line_to(LEFT_MARGIN + self.motif_end, y)
        context.stroke()

class Drawing:
    '''A drawing sequence.'''
    def __init__ (self, gene, exon, motif):
        self.gene = gene
        self.exon = exon
        self.motif = motif
        # self.context = context
    
    def draw(self, context):
        self.gene.draw_gene(context)
        self.gene.draw_exon(context)
        self.motif.draw_motif(context)

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

def motif_to_regex(motif:str) -> str:
    '''A funtion to return nucleotides correct nucletoides in motifs corresponding to the IUPAC dict.'''
    # eg motif= "ycgy", motif_regex = "[CT]cg[CT]" 
    # Using IUPAC dict
    iupac_dict = {
        "Y": "[CT]",
        "R": "[AG]",
        "S": "[GC]",
        "W": "[AT]",
        "K": "[GT]",
        "M": "[AC]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]",
        "N": "[ATCG]",
        "U": "[UT]",
    }

    # Replacing motifs with corresponding regex patterns
    
    regex_motif = motif.upper()
    for nucleotide, motif in iupac_dict.items():
        regex_motif = regex_motif.replace(nucleotide, motif)

    return regex_motif

# my_motif_regex = motif_to_regex("ygcy")
# print(my_motif_regex)

def motif_builder(gene_number: int, sequence: str, motifs: list[str]) -> list[Motif]:
    motif_list = []
    print(f"motif_builder {motifs=}")
    for motif_str in motifs:
        regex_motif = motif_to_regex(motif_str)
        # print(f'debug 888 {regex_motif=}')
        # motif_positions[motif] = [sequence, gene_number, for match in re.finditer(regex_motif, sequence)]
        for match in re.finditer(regex_motif, sequence, re.IGNORECASE):
            motif = Motif(match.start(), match.end(), gene_number, motif_str)
            motif_list.append(motif)
    return motif_list

def convert_motif_string_to_color(motif_str: str) -> tuple:
    '''Converts a motif string to RGB color'''
    # Define your color palette here
    color_palette = {
        "A": (0.9, 0.1, 0.1),  # Red
        "T": (0.1, 0.9, 0.1),  # Green
        "C": (0.1, 0.1, 0.9),  # Blue
        "G": (0.9, 0.9, 0.1)   # Yellow
        # Add more colors as needed
    }
    
    # Convert each nucleotide in the motif string to its corresponding color
    motif_color = [color_palette.get(nucleotide, (0, 0, 0)) for nucleotide in motif_str.upper()]
    
    # Average the colors to get a single color for the motif
    avg_color = tuple(sum(color[i] for color in motif_color) / len(motif_color) for i in range(3))
    
    return avg_color

# def draw_legend(context: cairo.Context, motif_colors: dict):
#     # Set font properties for legend labels
#     context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
#     context.set_font_size(12)

#     # Set initial position for legend
#     x = 650
#     y = 20

#     # Iterate over motif colors and draw legend
#     for motif, color in motif_colors.items():
#         # Draw colored rectangle
#         context.set_source_rgb(*color)
#         context.rectangle(x, y, 20, 10)
#         context.fill()

#         # Draw motif label
#         context.set_source_rgb(0, 0, 0)
#         context.move_to(x + 30, y + 10)
#         context.show_text(motif)

#         # Update position for next legend item
#         y += 20  # Increase y-coordinate for next legend item


# motif_list = motif_builder(1, "actgaaacccttgccttgggaaagcttgct", ["ygcy", "GCAUG"])
# print(motif_list)
# print("LJ is rad")

#############
# Arguments #
#############
            
def get_args():
    parser= argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Input fasta", required=True)
    parser.add_argument("-m", "--motifs", help="Input motifs", required=True)
    return parser.parse_args()

args = get_args()
fasta_filename=args.fasta
motif_filename=args.motifs

########
# Main #
########

motif_str_list = []
# # assign one color per motif
# # key = motifs (ygcy), values = colors (0.9, 0.1, 0.1)
# MOTIF_COLOR_DICT: dict[str, tuple[float, float, float]] = {} 
# COLOR_PALETTE = [
#         (0.9, 0.1, 0.1),  # Red
#         (0.1, 0.9, 0.1),  # Green
#         (0.1, 0.1, 0.9),  # Blue
#         (0.9, 0.9, 0.1),   # Yellow
#         (0.8, 0.3, 0.4),   # Mystery
#         # Add more colors as needed
# ]

# ygcy
# GCAUG
# catag
# YYYYYYYYYY

with open(motif_filename) as f:
    #print(f"debug 159 {motif_filename}")
    for index, line in enumerate(f):
        motif_str = line.strip()
        print(index)
        MOTIF_COLOR_DICT[motif_str] = COLOR_PALETTE[index]
        motif_str_list.append(motif_str)

oneline_fasta(fasta_filename, "intermediate.fasta")

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
    #print(f"debug 283 yes")
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
    #print(f"debug 352 {str(gene)=}")
    gene.draw_gene(context, gene_symbol)
    # gene_number += 1

    exon_start_end = exon_finder(sequence)
    exon_start = exon_start_end[0]
    exon_end = exon_start_end[1]

    # Making the exon object
    exon = Exon(exon_start, exon_end, gene_number)
    exon.draw_exon(context)

    motif_list = motif_builder(gene_number, sequence, motif_str_list)
    for motif in motif_list:
        # print(f"debug 111 {motif=}")
        motif.draw_motif(context)

# print(motif_list)
# print("LJ is rad")

surface.write_to_png("gene_sequence2.png")