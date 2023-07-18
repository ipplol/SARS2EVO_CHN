using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Net.Mail;

namespace CalLineageDistribution
{
    public class Month
    {
        public string monthName;
        public int totalSeqNumber = 0;
    }
    public class Lineage
    {
        public string lineageName;
        public Dictionary<string, int> Dic_lineageMonthCount = new Dictionary<string, int>();
    }
    class Program
    {
        static Dictionary<string, Month> Dic_month = new Dictionary<string, Month>();//存储名称和月份 列表
        static Dictionary<string, Lineage> Dic_lineage = new Dictionary<string, Lineage>();//储存名称和lineage 列表
        static void Main(string[] args)
        {
            StreamReader read = new StreamReader("M://China220701_230531/Data/allCHN.metadata.tsv");
            StreamReader readLineage = new StreamReader("M://China220701_230531/CHNSubtree/CHNLineage/LineageList.txt");
            StreamWriter write = new StreamWriter("M://China220701_230531/CHNSubtree/CHNLineage/LineageDistribution.txt");
            string line = readLineage.ReadLine();//读取要查找的lineage
            while(line!=null)
            {
                Lineage a = new Lineage();
                a.lineageName = line;
                Dic_lineage.Add(line, a);
                line = readLineage.ReadLine();
            }

            line = read.ReadLine();
            line = read.ReadLine();
            while (line!=null)
            {
                string[] line1 = line.Split('\t');
                if(line1[4].Length!=10)//日期不符合要求
                {
                    line = read.ReadLine();
                    continue;
                }

                /*if (!line1[1].Contains("EPI_ISL") && !line1[3].Contains("EPI_ISL"))//非GISAID序列
                {
                    line = read.ReadLine();
                    continue;
                }*/

                /*if (!line1[11].Contains("China"))//国家不符合要求
                {
                    line = read.ReadLine();
                    continue;
                }*/

                string date = line1[4].Substring(0, 7);

                //该月总序列数
                if(Dic_month.ContainsKey(date))//该月份存在，总数+1
                {
                    Dic_month[date].totalSeqNumber++;
                }
                else//添加新月份
                {
                    Month b = new Month();
                    b.monthName = date;
                    b.totalSeqNumber = 1;
                    Dic_month.Add(date, b);
                }

                //该lineage序列数
                if (Dic_lineage.ContainsKey(line1[18]))//该序列来自我们要找的lineage
                {
                    if (Dic_lineage[line1[18]].Dic_lineageMonthCount.ContainsKey(date))
                    { 
                        Dic_lineage[line1[18]].Dic_lineageMonthCount[date]++; 
                    }
                    else//该lineage还没有在这个月出现过
                    {
                        Dic_lineage[line1[18]].Dic_lineageMonthCount.Add(date, 1);
                    }
                }
                else
                {
                    string[] lineageFirst = line1[18].Split('.');
                    if (Dic_lineage.ContainsKey(lineageFirst[0]+"*"))//该序列来自我们要找的lineage 带*表示向下兼容子代 比如 A*
                    {
                        if (Dic_lineage[lineageFirst[0] + "*"].Dic_lineageMonthCount.ContainsKey(date))
                        {
                            Dic_lineage[lineageFirst[0] + "*"].Dic_lineageMonthCount[date]++;
                        }
                        else//该lineage还没有在这个月出现过
                        {
                            Dic_lineage[lineageFirst[0] + "*"].Dic_lineageMonthCount.Add(date, 1);
                        }
                    }
                }

                line = read.ReadLine();
            }

            string output = "Month\tSeqNum";
            foreach(var val in Dic_lineage)
            {
                output += "\t" + val.Key;
            }
            write.WriteLine(output);
            foreach (var val in Dic_month)
            {
                output = val.Key + "\t" + Convert.ToString(val.Value.totalSeqNumber);
                foreach(var valj in Dic_lineage)
                {
                    if(valj.Value.Dic_lineageMonthCount.ContainsKey(val.Key))
                    {
                        output += "\t" + Convert.ToString(Convert.ToDouble(valj.Value.Dic_lineageMonthCount[val.Key]) / val.Value.totalSeqNumber);
                    }
                    else
                    {
                        output += "\t0";
                    }
                }
                write.WriteLine(output);
            }

            readLineage.Close();
            read.Close();
            write.Close();
        }
    }
}
