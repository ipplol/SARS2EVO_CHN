using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MutincidenceGenomeDistribution
{
    public class Lineage
    {
        public string Name;
        public Dictionary<string, int> GeneNSMutCount = new Dictionary<string, int>();
        public Dictionary<string, int> GeneAllMutCount = new Dictionary<string, int>();
    }
    class Program
    {
        static void Main(string[] args)
        {
            List<Lineage> LineageList = new List<Lineage>();
            List<string> filelist = new List<string>();
            var files = Directory.GetFiles("M://China220701_230531/CHNSubtree/MutIncidence", "*.MutIncidence");
            foreach (string val in files)
                filelist.Add(val);
            int i, j, k;
            for(i=0;i<filelist.Count();i++)
            {
                StreamReader read = new StreamReader(filelist[i]);
                string[] line1 = filelist[i].Split('\\');
                string lineageName = line1[1].Substring(0, line1[line1.Length - 1].Length - 13);
                Lineage newl = new Lineage();
                newl.Name = lineageName;
                newl.GeneNSMutCount.Add("nsp1", 0);
                newl.GeneNSMutCount.Add("nsp2", 0);
                newl.GeneNSMutCount.Add("nsp3", 0);
                newl.GeneNSMutCount.Add("nsp4", 0);
                newl.GeneNSMutCount.Add("nsp5", 0);
                newl.GeneNSMutCount.Add("nsp6", 0);
                newl.GeneNSMutCount.Add("nsp7", 0);
                newl.GeneNSMutCount.Add("nsp8", 0);
                newl.GeneNSMutCount.Add("nsp9", 0);
                newl.GeneNSMutCount.Add("nsp10", 0);
                newl.GeneNSMutCount.Add("nsp11", 0);
                newl.GeneNSMutCount.Add("nsp12", 0);
                newl.GeneNSMutCount.Add("nsp13", 0);
                newl.GeneNSMutCount.Add("nsp14", 0);
                newl.GeneNSMutCount.Add("nsp15", 0);
                newl.GeneNSMutCount.Add("nsp16", 0);
                //newl.GeneNSMutCount.Add("ORF1ab", 0);
                newl.GeneNSMutCount.Add("ORF3a", 0);
                newl.GeneNSMutCount.Add("ORF6", 0);
                newl.GeneNSMutCount.Add("ORF7a", 0);
                newl.GeneNSMutCount.Add("ORF7b", 0);
                newl.GeneNSMutCount.Add("ORF8", 0);
                newl.GeneNSMutCount.Add("ORF10", 0);
                newl.GeneNSMutCount.Add("E", 0);
                newl.GeneNSMutCount.Add("M", 0);
                newl.GeneNSMutCount.Add("N", 0);
                newl.GeneNSMutCount.Add("S", 0);
                newl.GeneNSMutCount.Add("S_RBD", 0);
                newl.GeneNSMutCount.Add("S_NTD", 0);
                newl.GeneNSMutCount.Add("S_Other", 0);

                newl.GeneAllMutCount.Add("nsp1", 0);
                newl.GeneAllMutCount.Add("nsp2", 0);
                newl.GeneAllMutCount.Add("nsp3", 0);
                newl.GeneAllMutCount.Add("nsp4", 0);
                newl.GeneAllMutCount.Add("nsp5", 0);
                newl.GeneAllMutCount.Add("nsp6", 0);
                newl.GeneAllMutCount.Add("nsp7", 0);
                newl.GeneAllMutCount.Add("nsp8", 0);
                newl.GeneAllMutCount.Add("nsp9", 0);
                newl.GeneAllMutCount.Add("nsp10", 0);
                newl.GeneAllMutCount.Add("nsp11", 0);
                newl.GeneAllMutCount.Add("nsp12", 0);
                newl.GeneAllMutCount.Add("nsp13", 0);
                newl.GeneAllMutCount.Add("nsp14", 0);
                newl.GeneAllMutCount.Add("nsp15", 0);
                newl.GeneAllMutCount.Add("nsp16", 0);
                //newl.GeneAllMutCount.Add("ORF1ab", 0);
                newl.GeneAllMutCount.Add("ORF3a", 0);
                newl.GeneAllMutCount.Add("ORF6", 0);
                newl.GeneAllMutCount.Add("ORF7a", 0);
                newl.GeneAllMutCount.Add("ORF7b", 0);
                newl.GeneAllMutCount.Add("ORF8", 0);
                newl.GeneAllMutCount.Add("ORF10", 0);
                newl.GeneAllMutCount.Add("E", 0);
                newl.GeneAllMutCount.Add("M", 0);
                newl.GeneAllMutCount.Add("N", 0);
                newl.GeneAllMutCount.Add("S", 0);
                newl.GeneAllMutCount.Add("S_RBD", 0);
                newl.GeneAllMutCount.Add("S_NTD", 0);
                newl.GeneAllMutCount.Add("S_Other", 0);
                string line = read.ReadLine();
                line = read.ReadLine();
                while(line!=null)
                {
                    if(line.Contains("intergenic"))
                    {
                        line = read.ReadLine();
                        continue;
                    }
                    string[] line2 = line.Split('\t');
                    string[] line3 = line2[5].Split('_');
                    string gene = line3[0];
                    //if (gene.Contains("nsp")) gene = "ORF1ab";
                    if (line3[1][0] != line3[1][line3[1].Length - 1])//NS Mutation
                    {
                        newl.GeneNSMutCount[gene] += Convert.ToInt32(line2[6]);
                        if(gene == "S")
                        {
                            int pos = Convert.ToInt32(line3[1].Substring(1, line3[1].Length - 2));
                            if(pos >= 13 && pos <= 305)
                                newl.GeneNSMutCount["S_NTD"] += Convert.ToInt32(line2[6]);
                            else if (pos >= 331 && pos <= 531)
                                newl.GeneNSMutCount["S_RBD"] += Convert.ToInt32(line2[6]);
                            else
                                newl.GeneNSMutCount["S_Other"] += Convert.ToInt32(line2[6]);
                        }
                    }
                    
                        newl.GeneAllMutCount[gene] += Convert.ToInt32(line2[6]);
                        if (gene == "S")
                        {
                            int pos = Convert.ToInt32(line3[1].Substring(1, line3[1].Length - 2));
                            if (pos >= 13 && pos <= 305)
                                newl.GeneAllMutCount["S_NTD"] += Convert.ToInt32(line2[6]);
                            else if (pos >= 331 && pos <= 531)
                                newl.GeneAllMutCount["S_RBD"] += Convert.ToInt32(line2[6]);
                            else
                                newl.GeneAllMutCount["S_Other"] += Convert.ToInt32(line2[6]);
                        }
                    
                    line = read.ReadLine();
                }
                LineageList.Add(newl);
                read.Close();
            }

            StreamWriter write = new StreamWriter("M://China220701_230531/NSMutationGeneDistribution.txt");
            StreamWriter write2 = new StreamWriter("M://China220701_230531/AllMutationGeneDistribution.txt");
            //write.WriteLine("Lineage\tORF1ab\tSpike\tORF3a\tE\tM\tORF6\tORF7a\tORF7b\tORF8\tN\tORF10\tRBD\tNTD\tOther");
            //write2.WriteLine("Lineage\tORF1ab\tSpike\tORF3a\tE\tM\tORF6\tORF7a\tORF7b\tORF8\tN\tORF10\tRBD\tNTD\tOther");
            //string[] genes = "ORF1ab\tS\tORF3a\tE\tM\tORF6\tORF7a\tORF7b\tORF8\tN\tORF10\tS_RBD\tS_NTD\tS_Other".Split('\t');
            write.WriteLine("Lineage\tnsp1\tnsp2\tnsp3\tnsp4\tnsp5\tnsp6\tnsp7\tnsp8\tnsp9\tnsp10\tnsp11\tnsp12\tnsp13\tnsp14\tnsp15\tnsp16\tSpike\tORF3a\tE\tM\tORF6\tORF7a\tORF7b\tORF8\tN\tORF10\tRBD\tNTD\tOther");
            write2.WriteLine("Lineage\tnsp1\tnsp2\tnsp3\tnsp4\tnsp5\tnsp6\tnsp7\tnsp8\tnsp9\tnsp10\tnsp11\tnsp12\tnsp13\tnsp14\tnsp15\tnsp16\tSpike\tORF3a\tE\tM\tORF6\tORF7a\tORF7b\tORF8\tN\tORF10\tRBD\tNTD\tOther");
            string[] genes = "nsp1\tnsp2\tnsp3\tnsp4\tnsp5\tnsp6\tnsp7\tnsp8\tnsp9\tnsp10\tnsp11\tnsp12\tnsp13\tnsp14\tnsp15\tnsp16\tS\tORF3a\tE\tM\tORF6\tORF7a\tORF7b\tORF8\tN\tORF10\tS_RBD\tS_NTD\tS_Other".Split('\t');
            for (i=0;i<LineageList.Count();i++)
            {
                string output = LineageList[i].Name;
                for(j=0;j<genes.Length;j++)
                {
                    output += "\t" + LineageList[i].GeneNSMutCount[genes[j]];
                }
                write.WriteLine(output);
                output = LineageList[i].Name;
                for (j = 0; j < genes.Length; j++)
                {
                    output += "\t" + LineageList[i].GeneAllMutCount[genes[j]];
                }
                write2.WriteLine(output);
            }
            write.Close();
            write2.Close();
            return;
        }
    }
}
