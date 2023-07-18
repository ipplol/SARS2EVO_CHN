using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 从mutevent统计各lineage突变事件总数
{
    class Program
    {
        static void Main(string[] args)
        {
            Dictionary<string, int> EventCount = new Dictionary<string, int>();
            int i, j, k;
            StreamReader read = new StreamReader("//NAS8500/g/ChinaCoV/China220701_230531/Data/global_assignments.json.mutevent");
            StreamWriter write = new StreamWriter("//NAS8500/g/ChinaCoV/China220701_230531/Lineage_muteventCount.tsv");
            write.WriteLine("Lineage\tCount");
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if (EventCount.ContainsKey(line1[4]))
                    EventCount[line1[4]]++;
                else
                    EventCount.Add(line1[4], 1);
                line = read.ReadLine();
            }
            foreach (string val in EventCount.Keys)
                write.WriteLine(val + "\t" + EventCount[val]);
            read.Close();
            write.Close();
        }
    }
}
