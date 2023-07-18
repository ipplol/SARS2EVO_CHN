using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace SpikeOneStepMut
{
    class Program
    {
        static void Main(string[] args)
        {
            int i, j, k;
            string Sequence;
            Dictionary<string, string> Mimazi_Dic = new Dictionary<string, string>();
            StreamReader read = new StreamReader("M://China220701_230531/Group12/ACE2/BA5.Spike.fa");
            StreamWriter write = new StreamWriter("M://China220701_230531/Group12/ACE2/BA5.Spike.AvailableAA");
            StreamWriter writeNuc = new StreamWriter("M://China220701_230531/Group12/ACE2/BA5.Spike.AvailableNucAA");
            string line = read.ReadLine();
            Sequence = read.ReadLine();
            read.Close();

            read = new StreamReader("M://China220701_230531/Data/mimazi.txt");
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                Mimazi_Dic.Add(line1[0], line1[1]);
                line = read.ReadLine();
            }
            read.Close();

            for(i=0;i<Sequence.Length;i+=3)
            {
                string refCodon = "" + Sequence[i] + Sequence[i+1] + Sequence[i+2];
                string ATCG = "ATCG";
                List<string> mutAA = new List<string>();
                for(j=0;j<3;j++)
                {
                    for(k=0;k<4;k++)
                    {
                        char[] mutCodon = refCodon.ToCharArray();
                        mutCodon[j] = ATCG[k];
                        string mutCodon1 = new string(mutCodon);
                        if (Mimazi_Dic[mutCodon1] != Mimazi_Dic[refCodon])
                        {
                            write.WriteLine(Convert.ToString(i / 3 + 1) + Mimazi_Dic[mutCodon1]);
                            writeNuc.WriteLine(refCodon + "\t" + mutCodon1 + "\t" + refCodon[j] + "\t" + mutCodon1[j] + "\t" + Convert.ToString(i / 3 + 1) + Mimazi_Dic[mutCodon1]);
                        }
                    }
                }
            }
            write.Close();
            writeNuc.Close();
        }
    }
}
