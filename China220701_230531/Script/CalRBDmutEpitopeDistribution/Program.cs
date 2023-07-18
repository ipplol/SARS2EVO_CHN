using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace CalRBDmutEpitopeDistribution
{
    public class Mutation
    {
        public string MutName;
        public List<double> EscScores = new List<double>();
        public List<int> Esc12TF = new List<int>();//是否12类抗体的逃逸突变，1是，0否
        public int TotalEscEpitope = 0;//一共逃逸多少表位
        public bool NotEscMut = true;//是否为非逃逸突变
    }
    public class Clade
    {
        public string CladeName;
        public int CladeRBDTotalEvent = 0;
        public int TotalEscapeMutation = 0;
        public int TotalPositiveDMSMutation = 0;
        public int TotalESC_DMSMutation = 0;
        public List<double> Esc12Group = new List<double>();
    }
    class Program
    {
        static string workfold = "M://China220701_230531";
        static Dictionary<string, Mutation> MutationList = new Dictionary<string, Mutation>();
        static List<Clade> CladeList = new List<Clade>();
        static void WriteResult()
        {
            StreamWriter write = new StreamWriter(workfold + "/Group12/CladeRBDMut12Group.txt");
            StreamWriter writeCol = new StreamWriter(workfold + "/Group12/CladeRBDMut12Group.col.txt");
            StreamWriter writeTotalEvent = new StreamWriter(workfold + "/Group12/CladeRBDMutNumber.txt");
            int i, j, k;
            string[] epitopes = "Clade\tA1\tA2\tB\tC/D1\tD2\tD3\tD4\tE1/E2.1\tE2.2\tE3\tF1\tF3".Split('\t');
            write.WriteLine("Clade\tA1\tA2\tB\tC/D1\tD2\tD3\tD4\tE1/E2.1\tE2.2\tE3\tF1\tF3");
            writeCol.WriteLine("Clade\tEpitope\tEscProp");
            for (i = 0; i < CladeList.Count; i++)
            {
                string output = CladeList[i].CladeName;
                writeTotalEvent.WriteLine(output + "\t" + Convert.ToString(CladeList[i].CladeRBDTotalEvent));
                for (j = 0; j < 12; j++)
                {
                    output += "\t" + Convert.ToString((double)CladeList[i].Esc12Group[j] / CladeList[i].CladeRBDTotalEvent);
                    writeCol.WriteLine(CladeList[i].CladeName + "\t" + epitopes[j+1] + "\t" + Convert.ToString((double)CladeList[i].Esc12Group[j] / CladeList[i].CladeRBDTotalEvent));
                }
                write.WriteLine(output);
                
            }
            write.Close();
            writeCol.Close();
            writeTotalEvent.Close();

            write = new StreamWriter(workfold + "/Group12/CladeRBDEscMutProp.txt");
            write.WriteLine("Clade\tTotalEvent\tEscMut");
            for (i = 0; i < CladeList.Count; i++)
            {

                write.WriteLine(CladeList[i].CladeName + "\t" + Convert.ToString(CladeList[i].CladeRBDTotalEvent) + "\t" + Convert.ToString(CladeList[i].TotalEscapeMutation));
            }
            write.Close();

            return;
        }
        static void EscapeDistribution()
        {
            int i, j, k;
            List<string> filelist = new List<string>();
            var files = Directory.GetFiles(workfold + "/CHNSubtree/MutIncidence", "*.MutIncidence.RBDAA");
            foreach (string val in files)
                filelist.Add(val);

            for(i=0;i<filelist.Count;i++)
            {
                string[] line1 = filelist[i].Split('\\');

                Clade newc = new Clade();
                newc.CladeName = line1[1].Substring(0, line1[1].Length - 19);
                for (k = 0; k < 12; k++) newc.Esc12Group.Add(0);
                newc.CladeRBDTotalEvent = 0;
                newc.TotalEscapeMutation = 0;
                StreamReader read = new StreamReader(filelist[i]);
                string[] filename = filelist[i].Split('\\');
                StreamWriter write = new StreamWriter(workfold + "/Group12/ESC/" + filename[1] + ".A2");
                write.WriteLine("Clade\tMut\tIncidence\tA1\tA2\tB\tC/D1\tD2\tD3\tD4\tE1/E2.1\tE2.2\tE3\tF1\tF3");

                string line = read.ReadLine();
                line = read.ReadLine();
                while(line!=null)
                {
                    string[] line2 = line.Split('\t');
                    newc.CladeRBDTotalEvent += Convert.ToInt32(line2[1]);
                    if (MutationList.ContainsKey(line2[0].Substring(1, 4)) && !MutationList[line2[0].Substring(1, 4)].NotEscMut)//逃逸突变
                    {
                        newc.TotalEscapeMutation += Convert.ToInt32(line2[1]);
                        for (j = 0; j < 12; j++)
                        {
                            newc.Esc12Group[j] += MutationList[line2[0].Substring(1, 4)].Esc12TF[j] * Convert.ToInt32(line2[1]);
                            //newc.Esc12Group[j] += MutationList[line2[0].Substring(1, 4)].Esc12TF[j] * Convert.ToInt32(line2[1]) / MutationList[line2[0].Substring(1, 4)].TotalEscEpitope;
                        }
                        if(MutationList[line2[0].Substring(1, 4)].Esc12TF[1] > 0)//A2类逃逸
                        {
                            string output = filename[1] + "\t" + line2[0].Substring(1, 4) + "\t" + line2[1];
                            for (j = 0; j < 12; j++)
                            output += "\t" + MutationList[line2[0].Substring(1, 4)].Esc12TF[j];
                            write.WriteLine(output);
                        }
                    }
                    line = read.ReadLine();
                }
                read.Close();
                write.Close();
                CladeList.Add(newc);
            }
            return;
        }
        static void ReadinAndCalAverage3()
        {
            int mutationCount = 0;
            List<double> ScoreSum = new List<double>();
            StreamReader read = new StreamReader(workfold + "/Data/EscapeScore_PKU_NEW_BA.5.12EscOnly.txt"); //EscapeScore_PKU_AllWT.12EscOnly.txt
            StreamWriter write = new StreamWriter(workfold + "/Data/MutationEsc12_NEW_BA5.txt");
            //StreamWriter writed = new StreamWriter(workfold + "/EscapeMutationCal/EscScoreDistribution.txt");
            write.WriteLine("Mutation\tA1\tA2\tB\tC/D1\tD2\tD3\tD4\tE1/E2.1\tE2.2\tE3\tF1\tF3");
            //writed.WriteLine("Mutation\tEpitope\tEscapeScore\tType");
            string[] epitopelist = "A1\tA2\tB\tC/D1\tD2\tD3\tD4\tE1/E2.1\tE2.2\tE3\tF1\tF3".Split('\t');
            int i, j, k;
            for (i = 1; i <= 12; i++) ScoreSum.Add(0);
            string line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                Mutation newa = new Mutation();
                newa.MutName = line1[0];
                for (i = 1; i < line1.Length; i++)
                {
                    newa.EscScores.Add(Convert.ToDouble(line1[i]));
                    ScoreSum[i - 1] += Convert.ToDouble(line1[i]);
                }
                MutationList.Add(newa.MutName, newa);
                mutationCount++;
                line = read.ReadLine();
            }
            for (i = 0; i < 12; i++)
            {
                ScoreSum[i] = ScoreSum[i] / mutationCount;

                //阈值为3倍平均数
                ScoreSum[i] = ScoreSum[i]  * 3;

                //阈值为平均数+3倍方差
                /*double variance = 0;
                foreach (string val in MutationList.Keys)
                    variance += (MutationList[val].EscScores[i] - (ScoreSum[i])) * (MutationList[val].EscScores[i] - (ScoreSum[i]));
                variance /= mutationCount;
                ScoreSum[i] += 3 * variance;*/
            }
            

            foreach (string val in MutationList.Keys)
            {
                string output = val;
                for (j = 0; j < 12; j++)
                {
                    if (MutationList[val].EscScores[j] > ScoreSum[j])
                    {
                        MutationList[val].Esc12TF.Add(1);
                        MutationList[val].NotEscMut = false;
                        MutationList[val].TotalEscEpitope++;
                        output += "\t1";
                    }
                    else
                    {
                        MutationList[val].Esc12TF.Add(0);
                        output += "\t0";
                    }
                    //writed.WriteLine(val + "\t" + epitopelist[j] + "\t" + MutationList[val].EscScores[j] + "\tRBDMut");
                }
                write.WriteLine(output);
            }

            /*for (j = 0; j < 12; j++)
            {
                List<double> escscore = new List<double>();
                foreach (string val in MutationList.Keys)
                    escscore.Add(MutationList[val].EscScores[j]);
                escscore.Sort();
                for (i = 0; i < escscore.Count; i++)
                    writed.WriteLine(Convert.ToString(i) + "\t" + epitopelist[j] + "\t" + Convert.ToString(escscore[i]) + "\tRBDMut");
                writed.WriteLine(Convert.ToString(0) + "\t" + epitopelist[j] + "\t" + Convert.ToString(ScoreSum[j]) + "\tThreshold");
                writed.WriteLine(Convert.ToString(i) + "\t" + epitopelist[j] + "\t" + Convert.ToString(ScoreSum[j]) + "\tThreshold");
            }
            writed.Close();*/
            
            read.Close();
            write.Close();
            
            return;
        }
        static void Main(string[] args)
        {
            ReadinAndCalAverage3();//读入突变逃逸得分，计算突变是否逃逸12类抗体
            EscapeDistribution();//统计分布
            WriteResult();
            return;
        }
    }
}
