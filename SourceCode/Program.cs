using System;
using System.IO;
using System.Reflection;

namespace CSEP590A.HW5.ProteinCodingGenes
{
    /// <summary>
    /// Entry point for Program DiBoyed_ProteinCDS.exe
    /// Author      : Dipak Boyed
    /// Description : This program predicts protein coding sequence on a
    ///               a given DNA sequence. It outputs an ORF Length histogram 
    ///               of matching genes. The matches genes are predicted using
    ///               (a)stop codons locaion from a gene bank, and 
    ///               (b)a 3rd order Markov Model.
    /// </summary>
    class Program
    {
        #region static Member
        private static string FASTAFile = String.Empty;         // Input DNA sequence file 
        private static string GeneBankFile = String.Empty;      // Input GeneBank file
        #endregion

        #region static Methods
        /// <summary>
        ///  Method       : ShowUsage
        ///  Author       : Dipak Boyed
        ///  Description  : Prints Help menu
        /// </summary>
        static void ShowUsage()
        {
            Console.WriteLine();
            Console.WriteLine("Help Menu: '{0}'", Path.GetFileName(Assembly.GetExecutingAssembly().Location));
            Console.WriteLine();
            Console.WriteLine("This program predicts protein coding sequence on a given DNA sequence in FASTA format.");
            Console.WriteLine("Outputs an ORF Length histogram of matching genes. The gene matches are predicted using");
            Console.WriteLine(" (a) matching stop codons locations in a gene bank, and ");
            Console.WriteLine(" (b) a 3rd order Markov Model with training data for a trusted/background model coming");
            Console.WriteLine("     from ORFs of length greater/lesser than given thresholds.");
            Console.WriteLine("...");
            Console.WriteLine();
            Console.WriteLine("Usage: ");
            Console.WriteLine("{0} <sequence> <gene bank file>", Path.GetFileName(Assembly.GetExecutingAssembly().Location));
            Console.WriteLine();
            Console.WriteLine("Options: ");
            Console.WriteLine(" <sequence>          FASTA file representing the DNA sequence.");
            Console.WriteLine(" <gene bank file>]   GBK file for comparsion purposes only.");
            Console.WriteLine();
        }

        /// <summary>
        ///  Method       : ValidateArguments
        ///  Author       : Dipak Boyed
        ///  Description  : Validate the arguments count and values passed by user.
        /// </summary>
        /// <param name="args">string array representing the arguments</param>
        /// <returns>True if arguments are valid, false otherwise.</returns>
        static bool ValidateArguments(string[] args)
        {
            if (args.Length != 2)
            {
                return false;
            }
            else
            {
                Program.FASTAFile = args[0];
                Program.GeneBankFile = args[1];
                return true;
            }
        }

        /// <summary>
        ///  Method       : Main
        ///  Author       : Dipak Boyed
        ///  Description  : Entry point of the program DiBoyed_ProteinCDS
        /// </summary>
        /// <param name="args">string array representing the arguments passed by user</param>
        static void Main(string[] args)
        {
            if (!Program.ValidateArguments(args))
            {
                Program.ShowUsage();
                return;
            }
            try
            {
                Console.WriteLine("Predicting Protein coding sequence...");
                CDSAlgorithm cds = new CDSAlgorithm(Program.FASTAFile, Program.GeneBankFile);
                cds.Compute();
                cds.Print();
            }
            catch (Exception e)
            {
                Console.WriteLine();
                Console.WriteLine("Error: {0}", e.Message);
                Console.WriteLine(e.StackTrace);
                Console.WriteLine();
            }
        }
        #endregion
    }
}
