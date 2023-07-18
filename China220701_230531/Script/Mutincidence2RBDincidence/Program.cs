using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Mutincidence2RBDincidence //合并相同的氨基酸突变 排序并输出NS 突变
{
    class Lineage
    {

    }
    class Program
    {
        static void Main(string[] args)
        {
            string Workfold = "M://China220701_230531/CHNSubtree/MutIncidence";
            List<string> filelist = new List<string>();
            var files = Directory.GetFiles(Workfold, "*.MutIncidence");
            foreach (string val in files)
                filelist.Add(val);
            int i, j, k;
            for (i = 0; i < filelist.Count(); i++)
            {
                StreamReader read = new StreamReader(filelist[i]);
                List<string> Mutlist = new List<string>();
                List<int> Incidencelist = new List<int>();
                string line = read.ReadLine();
                line = read.ReadLine();
                string[] linez = line.Split('\t');
                while(line!=null)
                {
                    string[] line1 = line.Split('\t');
                    if(line1[5].Contains("S_"))//Spike
                    {
                        string[] line2 = line1[5].Split('_');
                        int pos = Convert.ToInt32(line2[1].Substring(1, line2[1].Length - 2));
                        if(pos>=331 && pos<=531)//RBD
                        {
                            if(line2[1][0]!=line2[1][line2[1].Length - 1])//NS
                            {
                                if(Mutlist.Contains(line2[1]))
                                {
                                    Incidencelist[Mutlist.IndexOf(line2[1])] += Convert.ToInt32(line1[6]);
                                }
                                else
                                {
                                    Mutlist.Add(line2[1]);
                                    Incidencelist.Add(Convert.ToInt32(line1[6]));
                                }
                            }
                        }
                    }
                    line = read.ReadLine();
                }
                read.Close();
                
                //sort and output
                int tmpi;
                string tmps;
                for(j=0;j<Mutlist.Count;j++)
                {
                    for(k=j+1;k<Mutlist.Count;k++)
                    {
                        if(Incidencelist[j]<Incidencelist[k])
                        {
                            tmps = Mutlist[j];Mutlist[j] = Mutlist[k];Mutlist[k] = tmps;
                            tmpi = Incidencelist[j];Incidencelist[j] = Incidencelist[k];Incidencelist[k] = tmpi;
                        }
                    }
                }
                StreamWriter write = new StreamWriter(filelist[i] + ".RBDAA");
                write.WriteLine("Mut\t" + linez[0]);
                for(j = 0; j < Mutlist.Count; j++)
                {
                    write.WriteLine(Mutlist[j] + "\t" + Convert.ToInt32(Incidencelist[j]));
                }
                write.Close();
            }
        }
    }
}
