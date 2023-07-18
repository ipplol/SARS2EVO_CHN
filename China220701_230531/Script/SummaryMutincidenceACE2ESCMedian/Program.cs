using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace SummaryMutincidenceACE2ESCMedian
{
    class Program
    {
        public static Dictionary<string, double> BindExpr = new Dictionary<string, double>();
        public static Dictionary<string, double> EscScore = new Dictionary<string, double>();
        static void Main(string[] args)
        {
            int i, j, k;
            StreamReader readf = new StreamReader("M://China220701_230531/Data/bindExpr.txt");
            string lineL = readf.ReadLine();
            while (lineL != null)
            {
                string[] lineL1 = lineL.Split('\t');
                BindExpr.Add(lineL1[0], Convert.ToDouble(lineL1[1]));
                lineL = readf.ReadLine();
            }
            readf.Close();

            readf = new StreamReader("M://China220701_230531/Data/EscapeScore_PKU_NEW_BA.5.single.txt");
            lineL = readf.ReadLine();
            lineL = readf.ReadLine();
            while (lineL != null)
            {
                string[] lineL1 = lineL.Split('\t');
                EscScore.Add(lineL1[0], Convert.ToDouble(lineL1[1]));
                lineL = readf.ReadLine();
            }
            readf.Close();

            StreamWriter write = new StreamWriter("M://China220701_230531/ChinaVSAbroad/BA5SubLineage/ACE2ESCMedian.txt");
            write.WriteLine("Lineage\tMedianACE2\tP25ACE2\tP75ACE2\tMedianESC\tP25ESC\tP75ESC");
            string Workfold = "M://China220701_230531/CHNSubtree/MutIncidence";
            List<string> filelist = new List<string>();
            var files = Directory.GetFiles(Workfold, "*.MutIncidence.RBDAA");
            foreach (string val in files)
                filelist.Add(val);

            for (k = 0; k < filelist.Count(); k++)
            {
                StreamReader read = new StreamReader(filelist[k]);
                List<double> ACEList = new List<double>();
                List<double> ESCList = new List<double>();
                string line = read.ReadLine();
                string[] title = line.Split('\t');
                line = read.ReadLine();
                while(line!=null)
                {
                    string[] line1 = line.Split('\t');
                    if (BindExpr.ContainsKey(line1[0].Substring(1, 4)))
                    {
                        for (i = 0; i < Convert.ToInt32(line1[1]); i++)
                        {
                            ACEList.Add(BindExpr[line1[0].Substring(1, 4)]);
                        }
                    }
                    if (EscScore.ContainsKey(line1[0].Substring(1, 4)))
                    {
                        for (i = 0; i < Convert.ToInt32(line1[1]); i++)
                        {
                            ESCList.Add(EscScore[line1[0].Substring(1, 4)]);
                        }
                    }
                    line = read.ReadLine();
                }
                ACEList.Sort();
                ESCList.Sort();
                string output = title[1];
                output += "\t" + Convert.ToString(ACEList[ACEList.Count / 2]);
                output += "\t" + Convert.ToString(ACEList[ACEList.Count / 4]);
                output += "\t" + Convert.ToString(ACEList[ACEList.Count / 4 * 3]);

                output += "\t" + Convert.ToString(ESCList[ESCList.Count / 2]);
                output += "\t" + Convert.ToString(ESCList[ESCList.Count / 4]);
                output += "\t" + Convert.ToString(ESCList[ESCList.Count / 4 * 3]);
                write.WriteLine(output);
                read.Close();
            }
            write.Close();
        }
    }
}
