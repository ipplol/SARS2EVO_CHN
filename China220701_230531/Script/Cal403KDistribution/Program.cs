using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Cal403KDistribution
{
    public class Lineage
    {
        public string Name;
        public List<int> CollectionDate = new List<int>();
        public int SeqWith403 = 0;
        public int TotalSeq = 0;
        public int MutEventWith403 = 0;
        public int TotalEvent = 0;
    }
    class Program
    {
        public static Dictionary<string, Lineage> LineageDic = new Dictionary<string, Lineage>();
        static void Main(string[] args)
        {
            int i, j, k;
            StreamReader read = new StreamReader("M://China220701_230531/Data/global_assignments.json.nodeinfo");
            string line = read.ReadLine();
            line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                if (line.Substring(0, 5) != "node_")
                {
                    string[] line1 = line.Split('\t');
                    if (!LineageDic.ContainsKey(line1[3]))
                    {
                        Lineage newl = new Lineage();
                        newl.Name = line1[3];
                        LineageDic.Add(line1[3], newl);
                    }
                    LineageDic[line1[3]].TotalSeq++;
                    if (line1[4].Contains("G22770A") || line1[4].Contains("T22770A") || line1[4].Contains("C22770A"))
                        LineageDic[line1[3]].SeqWith403++;
                    if (line1[1].Length == 10)
                        LineageDic[line1[3]].CollectionDate.Add(Convert.ToInt32(line1[1].Substring(0, 4) + line1[1].Substring(5, 2) + line1[1].Substring(8, 2)));
                }
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader("M://China220701_230531/Data/global_assignments.json.mutevent");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                if (line.Substring(0, 5) != "node_")
                {
                    string[] line1 = line.Split('\t');
                    LineageDic[line1[4]].TotalEvent++;
                    if (line1[3].Contains("G22770A") || line1[3].Contains("T22770A") || line1[3].Contains("C22770A"))
                        LineageDic[line1[4]].MutEventWith403++;
                }
                line = read.ReadLine();
            }
            read.Close();

            StreamWriter write = new StreamWriter("M://China220701_230531/ChinaVSAbroad/403KDistribution/Global_History_AllLineage_403K_Distribution.tsv");
            write.WriteLine("Lineage\tCollectionDate5P\tCollectionDate5PYear\tTotalSeq\tThe403KSeq\tSeqProp\tTotalEvent\tThe403KEvent\tEventProp");
            foreach(string val in LineageDic.Keys)
            {
                if (val != "" && LineageDic[val].TotalSeq >= 20)
                {
                    string output = LineageDic[val].Name + "\t";
                    LineageDic[val].CollectionDate.Sort();
                    string date = "NA";
                    if(LineageDic[val].CollectionDate.Count > 0)
                        date = Convert.ToString(LineageDic[val].CollectionDate[Convert.ToInt32(LineageDic[val].CollectionDate.Count() * 0.05)]);
                    output += date + "\t";
                    if (date != "NA")
                        output += Convert.ToString(Convert.ToDouble(date.Substring(0, 4)) + Math.Min(1, ((Convert.ToDouble(date.Substring(4, 2)) - 1) * 30 + Convert.ToDouble(date.Substring(6, 2))) / 365)) + "\t";
                    else
                    {
                        output += "NA\t";
                        continue;
                    }
                    output += Convert.ToString(LineageDic[val].TotalSeq) + "\t";
                    output += Convert.ToString(LineageDic[val].SeqWith403) + "\t";
                    output += Convert.ToString(Convert.ToDouble(LineageDic[val].SeqWith403) / Convert.ToDouble(LineageDic[val].TotalSeq)) + "\t";

                    output += Convert.ToString(LineageDic[val].TotalEvent) + "\t";
                    output += Convert.ToString(LineageDic[val].MutEventWith403) + "\t";
                    output += Convert.ToString(Convert.ToDouble(LineageDic[val].MutEventWith403) / Convert.ToDouble(LineageDic[val].TotalEvent)) + "\t";

                    write.WriteLine(output);
                }
            }

            write.Close();
        }
    }
}
