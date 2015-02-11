using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace CSEP590A.HW5.ProteinCodingGenes
{
    /// <summary>
    /// Nucleotide Enumerator
    /// </summary>
    internal enum DNABase : int
    {
        A = 0,
        C = 1,
        G = 2,
        T = 3,
    }

    internal class ORF
    {
        public int Start;
        public int End;
        public bool isRealGene;
        private bool inReverseOrder;
        public double LogPQScore;

        internal ORF(int start, int end, bool isReal, bool inReverse)
        {
            this.Start = start;
            this.End = end;
            this.isRealGene = isReal;
            this.inReverseOrder = inReverse;
        }

        internal ORF(int start, int end, bool isReal) : this(start, end, isReal, false) {}

        public int Length
        {
            get
            {
                return (this.End - this.Start + 1);
            }
        }
    }

    /// <summary>
    ///  Class        : CDSAlgorithm
    ///  Author       : Dipak Boyed
    ///  Description  : class representing the Protein CDS Prediction algorithm
    /// </summary>
    internal class CDSAlgorithm
    {
        #region static Members
        public static int PThresholdLength = 1400;  // Threshold length for trusted markov model
        public static int QThresholdLength = 50;    // Threshold length for background markov model

        // Stop Codon sequence
        private static DNABase[,] StopCodons = new DNABase[,] {
        /* TAA*/ { DNABase.T, DNABase.A, DNABase.A },
        /* TAG*/ { DNABase.T, DNABase.A, DNABase.G },
        /* TGA*/ { DNABase.T, DNABase.G, DNABase.A }
        };
        #endregion

        #region private Members
        private List<DNABase> DNASequence;
        private Dictionary<int, List<ORF>> CodingSequence;

        private Dictionary<int, bool> StopCodonsFromGeneBank;

        private double[]       P1Values;         // P1[x1] == P(Nucleotide in position x1)
        private double[]       Q1Values;
        private double[,]      P2Values;         // P2[x1, x2] == P(x2|x1)
        private double[,]      Q2Values;
        private double[, ,]    P3Values;         // P3[x1, x2, x3] = P(x3|x1x2)
        private double[, ,]    Q3Values;
        private double[, , ,]  P4Values;         // P4[x1, x2, x3, x4] = P(x4|x1x2x3)
        private double[, , ,]  Q4Values;
        #endregion

        #region Constructors
        /// <summary>
        ///  Method       : Ctor
        ///  Author       : Dipak Boyed
        ///  Description  : Constructs an CDSAlgorithm instance.
        ///                 Also parses the DNA sequence from a FASTA file
        /// </summary>
        /// <param name="fastaFile">FASTA file name containing DNA sequence.</param>
        /// <param name="gbkFile"> GeneBank File for comparison</param>
        internal CDSAlgorithm(string fastaFile, string gbkFile)
        {
            this.DNASequence = new List<DNABase>();
            this.CodingSequence = new Dictionary<int, List<ORF>>();
            this.StopCodonsFromGeneBank = new Dictionary<int, bool>();

            this.ReadSequence(fastaFile);
            if (!String.IsNullOrEmpty(gbkFile))
            {
                this.ReadStopCodonsFromGeneBank(gbkFile);
            }

            this.P1Values = new double[4];
            this.P2Values = new double[4, 4];
            this.P3Values = new double[4, 4, 4];
            this.P4Values = new double[4, 4, 4, 4];
            this.Q1Values = new double[4];
            this.Q2Values = new double[4, 4];
            this.Q3Values = new double[4, 4, 4];
            this.Q4Values = new double[4, 4, 4, 4];
        }
        #endregion

        #region Helper Methods
        /// <summary>
        ///  Method       : ReadSequence
        ///  Author       : Dipak Boyed
        ///  Description  : Ensures a valid FASTA file is specified.
        ///                 Also ensures each character in the sequence is a valid DNA base.
        ///                 Reads DNA sequence into internal data structure.
        /// </summary>
        /// <param name="sequence">string representing FASTA file name.</param>
        private void ReadSequence(string fastaFile)
        {
            string validatedSequence = String.Empty;
            if (String.IsNullOrEmpty(fastaFile))
            {
                throw new ArgumentException("FASTA file name cannot be null or empty string.", fastaFile);
            }

            if (File.Exists(fastaFile))
            {
                // sequence specified as a FASTA file
                Console.WriteLine("     File '{0}' exists. Attempting to read FASTA file...", Path.GetFileName(fastaFile));
                int lineNumber = 0;
                Console.Write("Reading line number: {0,5}", lineNumber);
                using (StreamReader reader = new StreamReader(fastaFile))
                {
                    while (!reader.EndOfStream)
                    {
                        string line = reader.ReadLine();
                        if (Console.CursorLeft >= 5)
                        {
                            Console.SetCursorPosition(Console.CursorLeft - 5, Console.CursorTop);
                            Console.Write("{0,5}", ++lineNumber);
                        }
                        if (line.StartsWith(">"))
                        {
                            if (Console.CursorLeft >= 26)
                            {
                                Console.SetCursorPosition(Console.CursorLeft - 26, Console.CursorTop);
                                Console.WriteLine("     Ignoring comments:\t'{0}'", line);
                                Console.Write("Reading line number: {0,5}", lineNumber);
                            }
                        }
                        else
                        {
                            validatedSequence += line;
                        }
                    }
                }
                Console.WriteLine();
            }
            else
            {
               throw new ArgumentException(String.Format("FASTA file '{0}' not found. Must specify a valid FASTA file", Path.GetFullPath(fastaFile)));
            }

            // Ensure each character in the sequence is a valid amino acid
            for (int i = 0; i < validatedSequence.Length; i++)
            {
                DNABase currentBase = DNABase.T;
                switch(validatedSequence[i])
                {
                    case 'A':
                        currentBase = DNABase.A;
                        break;
                    case 'C':
                        currentBase = DNABase.C;
                        break;
                    case 'G':
                        currentBase = DNABase.G;
                        break;
                    case 'T':
                        currentBase = DNABase.T;
                        break;
                    default:
                        Console.WriteLine("Found unknown base '{0}' at position '{1}'. Treating it as base 'T'.", validatedSequence[i], i);
                        currentBase = DNABase.T;
                        break;
                }
                this.DNASequence.Add(currentBase);
            }
            Console.WriteLine("     Successfully read sequence of length {0}.", this.DNASequence.Count);
        }

        /// <summary>
        ///  Method       : ReadStopCodonsFromGeneBank
        ///  Author       : Dipak Boyed
        ///  Description  : Reads stop codon locations of re al genes specified
        ///                 in the gene bank.
        /// </summary>
        /// <param name="gbkFile">Gene bank file to parse</param>
        private void ReadStopCodonsFromGeneBank(string gbkFile)
        {
            if (String.IsNullOrEmpty(gbkFile))
            {
                throw new ArgumentException("Emission file name cannot be null or empty string.", gbkFile);
            }

            if (File.Exists(gbkFile))
            {
                Console.WriteLine("     Attempting to read Gene Bank file {0}. Stop codons in gene bank will be used for comparison...", Path.GetFileName(gbkFile));
                Regex CDSRegex = new Regex(@"CDS\s*\d+\.\.(\d+)", RegexOptions.Compiled);
                using (StreamReader reader = new StreamReader(gbkFile))
                {
                    while (!reader.EndOfStream)
                    {
                        string line = reader.ReadLine().Trim();
                        if (!String.IsNullOrEmpty(line))
                        {
                            if (CDSRegex.IsMatch(line))
                            {
                                this.StopCodonsFromGeneBank.Add(Convert.ToInt32(CDSRegex.Match(line).Groups[1].Value), true);
                            }
                        }
                    }
                }
            }
            else
            {
                throw new ArgumentException(String.Format("Emission file '{0}' not found. Must specify a valid file name.", Path.GetFullPath(gbkFile)));
            }
        }

        /// <summary>
        ///  Method       : IsStopCodon
        ///  Author       : Dipak Boyed
        ///  Description  : Determine if a given codon is a stop codon or not
        /// </summary>
        /// <param name="codon">DNA base sequeunce of length 3 (codon)</param>
        /// <returns>True if a stop codon, false otherwise</returns>
        private bool IsStopCodon(DNABase[] codon)
        {
            return ((CDSAlgorithm.StopCodons[0, 0] == codon[0] && CDSAlgorithm.StopCodons[0, 1] == codon[1] && CDSAlgorithm.StopCodons[0, 2] == codon[2]) ||
                    (CDSAlgorithm.StopCodons[1, 0] == codon[0] && CDSAlgorithm.StopCodons[1, 1] == codon[1] && CDSAlgorithm.StopCodons[1, 2] == codon[2]) ||
                    (CDSAlgorithm.StopCodons[2, 0] == codon[0] && CDSAlgorithm.StopCodons[2, 1] == codon[1] && CDSAlgorithm.StopCodons[2, 2] == codon[2]));
        }

        /// <summary>
        ///  Method       : ScanORFs
        ///  Author       : Dipak Boyed
        ///  Description  : Scan the DNA sequence in reading frames to find ORFs ending with a stop codon
        /// </summary>
        private void ScanORFs()
        {
            int i = 0;
            int orf1Start = i;
            int orf2Start = i + 1;
            int orf3Start = i + 2;
            for (; i < this.DNASequence.Count; i = i + 3)
            {
                // read frame 1 (which started at offset 0)
                if (i + 2 < this.DNASequence.Count)
                {
                    if (this.IsStopCodon(this.DNASequence.GetRange(i, 3).ToArray()))
                    {
                        ORF orf = new ORF(orf1Start + 1, i + 3, this.StopCodonsFromGeneBank.ContainsKey(i + 3) && this.StopCodonsFromGeneBank[i + 3]);
                        if (!this.CodingSequence.ContainsKey(orf.Length))
                        {
                            this.CodingSequence.Add(orf.Length, new List<ORF>());
                        }
                        this.CodingSequence[orf.Length].Add(orf);
                        orf1Start = i + 3;
                        // See if this CDS needs to be added to Markov Model training data
                        if (orf.Length > CDSAlgorithm.PThresholdLength)
                        {
                            this.AddToPTrainingData(orf);
                        }
                        else if (orf.Length < CDSAlgorithm.QThresholdLength)
                        {
                            this.AddToQTrainingData(orf);
                        }
                    }
                }

                // read frame 2 (which started at offset 1)
                if (i + 3 < this.DNASequence.Count)
                {
                    if (this.IsStopCodon(this.DNASequence.GetRange(i + 1, 3).ToArray()))
                    {
                        ORF orf = new ORF(orf2Start + 1, i + 4, this.StopCodonsFromGeneBank.ContainsKey(i + 4) && this.StopCodonsFromGeneBank[i + 4]);
                        if (!this.CodingSequence.ContainsKey(orf.Length))
                        {
                            this.CodingSequence.Add(orf.Length, new List<ORF>());
                        }
                        this.CodingSequence[orf.Length].Add(orf);
                        orf2Start = i + 4;
                        // See if this CDS needs to be added to Markov Model training data
                        if (orf.Length > CDSAlgorithm.PThresholdLength)
                        {
                            this.AddToPTrainingData(orf);
                        }
                        else if (orf.Length < CDSAlgorithm.QThresholdLength)
                        {
                            this.AddToQTrainingData(orf);
                        }
                    }
                }

                // read frame 3 (which started at offset 2)
                if (i + 4 < this.DNASequence.Count)
                {
                    if (this.IsStopCodon(this.DNASequence.GetRange(i + 2, 3).ToArray()))
                    {
                        ORF orf = new ORF(orf3Start + 1, i + 5, this.StopCodonsFromGeneBank.ContainsKey(i + 5) && this.StopCodonsFromGeneBank[i + 5]);
                        if (!this.CodingSequence.ContainsKey(orf.Length))
                        {
                            this.CodingSequence.Add(orf.Length, new List<ORF>());
                        }
                        this.CodingSequence[orf.Length].Add(orf);
                        orf3Start = i + 5;
                        // See if this CDS needs to be added to Markov Model training data
                        if (orf.Length > CDSAlgorithm.PThresholdLength)
                        {
                            this.AddToPTrainingData(orf);
                        }
                        else if (orf.Length < CDSAlgorithm.QThresholdLength)
                        {
                            this.AddToQTrainingData(orf);
                        }
                    }
                }
            }
        }

        private void AddToPTrainingData(ORF orf)
        {
            // subtract 1 from initial/end value as the DNA sequence is stored in a 0-index based list
            // subtract 3 from end value in order to EXCLUDE stop codon from training data
            for (int index = orf.Start - 1; index < orf.Start + orf.Length - 4; index++)
            {
                this.P1Values[(int)this.DNASequence[index]]++;
                if (index > orf.Start - 1)
                {
                    this.P2Values[(int)this.DNASequence[index - 1], (int)this.DNASequence[index]]++;
                    if (index > orf.Start)
                    {
                        this.P3Values[(int)this.DNASequence[index - 2], (int)this.DNASequence[index - 1], (int)this.DNASequence[index]]++;
                        if (index > orf.Start + 1)
                        {
                            this.P4Values[(int)this.DNASequence[index - 3],
                                          (int)this.DNASequence[index - 2],
                                          (int)this.DNASequence[index - 1],
                                          (int)this.DNASequence[index]]++;
                        }
                    }
                }
            }
        }

        private void AddToQTrainingData(ORF orf)
        {
            // subtract 1 from initial value as the DNA sequence is stored in a 0-index based list
            // subtract 3 from end value in order to EXCLUDE stop codon from training data
            for (int index = orf.Start - 1; index < orf.Start + orf.Length - 4; index++)
            {
                this.Q1Values[(int)this.DNASequence[index]]++;
                if (index > orf.Start - 1)
                {
                    this.Q2Values[(int)this.DNASequence[index - 1], (int)this.DNASequence[index]]++;
                    if (index > orf.Start)
                    {
                        this.Q3Values[(int)this.DNASequence[index - 2], (int)this.DNASequence[index - 1], (int)this.DNASequence[index]]++;
                        if (index > orf.Start + 1)
                        {
                            this.Q4Values[(int)this.DNASequence[index - 3],
                                          (int)this.DNASequence[index - 2],
                                          (int)this.DNASequence[index - 1],
                                          (int)this.DNASequence[index]]++;
                        }
                    }
                }
            }
        }

        private void ComputePQTrainingData()
        {
            // - Instance of training data were added to during the scan of ORFs
            // - Now divide those instances with total occurences to get probabilities
            // - Store log probabilities
            double totalP1Occurences = this.P1Values[0] + this.P1Values[1] + this.P1Values[2] + this.P1Values[3];
            double totalQ1Occurences = this.Q1Values[0] + this.Q1Values[1] + this.Q1Values[2] + this.Q1Values[3];
            for (int i = 0; i < 4; i++)
            {
                if (totalP1Occurences != 0) this.P1Values[i] = Math.Log(this.P1Values[i] / totalP1Occurences);
                if (totalQ1Occurences != 0) this.Q1Values[i] = Math.Log(this.Q1Values[i] / totalQ1Occurences);

                double totalP2_i_Occurences = this.P2Values[i, 0] + this.P2Values[i, 1] + this.P2Values[i, 2] + this.P2Values[i, 3];
                double totalQ2_i_Occurences = this.Q2Values[i, 0] + this.Q2Values[i, 1] + this.Q2Values[i, 2] + this.Q2Values[i, 3];

                for (int j = 0; j < 4; j++)
                {
                    if (totalP2_i_Occurences != 0) this.P2Values[i, j] = Math.Log(this.P2Values[i, j] / totalP2_i_Occurences);
                    if (totalQ2_i_Occurences != 0) this.Q2Values[i, j] = Math.Log(this.Q2Values[i, j] / totalQ2_i_Occurences);

                    double totalP3_ij_Occurences = this.P3Values[i, j, 0] + this.P3Values[i, j, 1] + this.P3Values[i, j, 2] + this.P3Values[i, j, 3];
                    double totalQ3_ij_Occurences = this.Q3Values[i, j, 0] + this.Q3Values[i, j, 1] + this.Q3Values[i, j, 2] + this.Q3Values[i, j, 3];

                    for (int k = 0; k < 4; k++)
                    {
                        if (totalP3_ij_Occurences != 0) this.P3Values[i, j, k] = Math.Log(this.P3Values[i, j, k] / totalP3_ij_Occurences);
                        if (totalQ3_ij_Occurences != 0) this.Q3Values[i, j, k] = Math.Log(this.Q3Values[i, j, k] / totalQ3_ij_Occurences);

                        double totalP4_ijk_Occurences = this.P4Values[i, j, k, 0] + this.P4Values[i, j, k, 1] + this.P4Values[i, j, k, 2] + this.P4Values[i, j, k, 3];
                        double totalQ4_ijk_Occurences = this.Q4Values[i, j, k, 0] + this.Q4Values[i, j, k, 1] + this.Q4Values[i, j, k, 2] + this.Q4Values[i, j, k, 3];

                        for (int l = 0; l < 4; l++)
                        {
                            if (totalP4_ijk_Occurences != 0) this.P4Values[i, j, k, l] = Math.Log(this.P4Values[i, j, k, l] / totalP4_ijk_Occurences);
                            if (totalQ4_ijk_Occurences != 0) this.Q4Values[i, j, k, l] = Math.Log(this.Q4Values[i, j, k, l] / totalQ4_ijk_Occurences);
                        }
                    }
                }
            }
        }

        private void ScoreORFs()
        {
            foreach (int key in this.CodingSequence.Keys)
            {
                for (int i = 0; i < this.CodingSequence[key].Count; i++)
                {
                    this.CodingSequence[key][i].LogPQScore = this.ComputeORFScore(this.CodingSequence[key][i]);
                }
            }
        }

        private double ComputeORFScore(ORF orf)
        {
            double LogPScore = 0.0;
            double LogQScore = 0.0;

            // subtract 1 from initial/end value as the DNA sequence is stored in a 0-index based list
            // subtract 3 from end value in order to EXCLUDE stop codon from training data
            for (int index = orf.Start - 1; index < orf.Start + orf.Length - 4; index++)
            {
                if (index == orf.Start - 1)
                {
                    LogPScore += this.P1Values[(int)this.DNASequence[index]];
                    LogQScore += this.Q1Values[(int)this.DNASequence[index]];
                }
                else if (index == orf.Start)
                {
                    LogPScore += this.P2Values[(int)this.DNASequence[index - 1], (int)this.DNASequence[index]];
                    LogQScore += this.Q2Values[(int)this.DNASequence[index - 1], (int)this.DNASequence[index]];
                }
                else if (index == orf.Start + 1)
                {
                    LogPScore += this.P3Values[(int)this.DNASequence[index - 2], (int)this.DNASequence[index - 1], (int)this.DNASequence[index]];
                    LogQScore += this.Q3Values[(int)this.DNASequence[index - 2], (int)this.DNASequence[index - 1], (int)this.DNASequence[index]];
                }
                else
                {
                    LogPScore += this.P4Values[(int)this.DNASequence[index - 3], (int)this.DNASequence[index - 2], (int)this.DNASequence[index - 1], (int)this.DNASequence[index]];
                    LogQScore += this.Q4Values[(int)this.DNASequence[index - 3], (int)this.DNASequence[index - 2], (int)this.DNASequence[index - 1], (int)this.DNASequence[index]];
                }
            }
            return LogPScore - LogQScore;
        }
        #endregion

        #region public Methods
        /// <summary>
        ///  Method       : Compute
        ///  Author       : Dipak Boyed
        ///  Description  : Compute the prediction of coding sequences
        /// </summary>
        public void Compute()
        {
            Console.WriteLine("    Scanning for ORFs...");
            // scan ORFs for coding sequence
            this.ScanORFs();

            Console.WriteLine("    Computing training data (P/Q) for Markov Model");
            Console.WriteLine("    Threshold Length for trusted model training data   : > {0}", CDSAlgorithm.PThresholdLength);
            Console.WriteLine("    Threshold Length for background model training data: < {0}", CDSAlgorithm.QThresholdLength);
            // for each ORF parse it's coding sequence to train P/Q distribution
            this.ComputePQTrainingData();

            Console.WriteLine("    Scoring [log(P(x)/Q(x)] each ORF based on training data");
            this.ScoreORFs();
        }

        /// <summary>
        ///  Method       : Print
        ///  Author       : Dipak Boyed
        ///  Description  : Print results of the coding sequence prediction.
        ///                 Prints an ORF Length Histogram containing:
        ///                   - [Count]  Count of  ORF of that given length
        ///                   - [#GBK]   No. of matching genes whose stop codon locations matched those specified in the GBK.
        ///                   - [#MM]    No. of matching genes predicted by the 3rd order trusted Markov Model
        ///                   - [#MMGBK] No. of matching #MM genes that also match #GBK.
        ///                   - [Avg.]   Average log(P(x)/Q(x)) of ORFs of the given length
        /// </summary>
        public void Print()
        {
            Console.WriteLine();
            StringBuilder stringBuilder = new StringBuilder(" Printing ORF Histogram...");
            stringBuilder.AppendLine();
            stringBuilder.AppendLine(" ----------------------------------------------------------------");
            stringBuilder.AppendLine(" NOTE: ORF Length is no. of nucleotides and includes stop codons");
            stringBuilder.AppendLine(" ----------------------------------------------------------------");
            stringBuilder.AppendLine(" Length Count  #GBK   #MM #MMGBK Avg.  ");
            stringBuilder.AppendLine(" ----------------------------------------------------------------");
            List<KeyValuePair<int, List<ORF>>> codingSequencesSortedByLength = new List<KeyValuePair<int,List<ORF>>>(this.CodingSequence);
            codingSequencesSortedByLength.Sort(
                delegate( KeyValuePair<int, List<ORF>>firstPair, KeyValuePair<int, List<ORF>> secondPair)
                {
                    return firstPair.Key.CompareTo(secondPair.Key);
                });
            int totalMatchFromGBK = 0;
            int totalNonMatchFromGBK = 0;
            int totalMatchFromMM = 0;
            int totalMatchFromBoth = 0;
            int totalFalseNegFromMM = 0;
            foreach (KeyValuePair<int, List<ORF>> keyValuePair in codingSequencesSortedByLength)
            {
                int gbkRealGenes = 0;
                int markovModelRealGenes = 0;
                int bothRealGenes = 0;
                double averageLogPQScore = 0.0;
                foreach (ORF orf in keyValuePair.Value)
                {
                    if (orf.isRealGene)
                        gbkRealGenes++;

                    if (orf.LogPQScore > 0)
                    {
                        markovModelRealGenes++;
                        if (orf.isRealGene)
                        {
                            bothRealGenes++;
                        }
                    }
                    else
                    {
                        if (orf.isRealGene)
                            totalFalseNegFromMM++;
                    }
                    averageLogPQScore += orf.LogPQScore;
                }
                averageLogPQScore = averageLogPQScore / (double)keyValuePair.Value.Count;
                string formattedString = "{0} : {1}  {2}  {3}  {4}  {5:0.00000}";
                if (keyValuePair.Key <= 9999 && keyValuePair.Value.Count <= 99999 && gbkRealGenes <= 99999 && markovModelRealGenes <= 99999 && bothRealGenes <= 99999)
                {
                    formattedString = " {0,4} : {1,5} {2,5} {3,5} {4,5}  {5:0.00000}";
                }
                stringBuilder.AppendFormat(formattedString, 
                                           keyValuePair.Key, 
                                           keyValuePair.Value.Count, 
                                           gbkRealGenes, 
                                           markovModelRealGenes, 
                                           bothRealGenes, 
                                           averageLogPQScore);
                stringBuilder.AppendLine();
                totalMatchFromGBK += gbkRealGenes;
                totalNonMatchFromGBK += (keyValuePair.Value.Count - gbkRealGenes);
                totalMatchFromMM += markovModelRealGenes;
                totalMatchFromBoth += bothRealGenes;
            }
            stringBuilder.AppendLine(" ----------------------------------------------------------------");
            stringBuilder.AppendLine(" #GBK Approach Summary:");
            stringBuilder.AppendLine(String.Format("    Total Genes specified in the GeneBank: {0}", this.StopCodonsFromGeneBank.Count));
            stringBuilder.AppendLine(String.Format("    Total ORFs found: {0}, Genes predicted in #GBK: {1}, Difference: {2}", totalMatchFromGBK + totalNonMatchFromGBK, totalMatchFromGBK, totalNonMatchFromGBK));
            stringBuilder.AppendLine();
            stringBuilder.AppendLine(" #MM Approach Summary:");
            stringBuilder.AppendFormat("    True Positives  ( #MM ^  #GBK): {0}", totalMatchFromBoth);
            stringBuilder.AppendLine();
            stringBuilder.AppendFormat("    False Positives ( #MM ^ !#GBK): {0}", totalMatchFromMM - totalMatchFromBoth);
            stringBuilder.AppendLine();
            stringBuilder.AppendFormat("    False Negatives (!#MM ^  #GBK): {0}", totalFalseNegFromMM);
            stringBuilder.AppendLine();
            stringBuilder.AppendLine(" ----------------------------------------------------------------");
            Console.Write(stringBuilder.ToString());
        }
        #endregion
    }
}
