# ProteinCodingGenes
http://courses.cs.washington.edu/courses/csep590a/13sp/hw/hw5.html 

HW5: Protein CDS Gene Prediction
---------------------------------------------------------------------

List of files submitted:
------------------------------
Report.doc/pdf  	: Report for Project5
SourceCode		: Source code directory
OutputLogs		: Directory containing output logs

What is this Project about?
------------------------------------
See the Report pdf for details

How to run this program?
-----------------------------------
F5 the Visual Studio solution to build and run

This program predicts protein coding sequence on a given DNA sequence in FASTA format.
Outputs an ORF Length histogram of matching genes. The gene matches are predicted using
 (a) matching stop codons locations in a gene bank, and
 (b) a 3rd order Markov Model with training data for a trusted/background model coming
     from ORFs of length greater/lesser than given thresholds.
...

Usage:
DiBoyed_ProteinCDS.exe <sequence> <gene bank file>

Options:
 <sequence>          FASTA file representing the DNA sequence.
 <gene bank file>    GBK file for comparsion purposes only.
