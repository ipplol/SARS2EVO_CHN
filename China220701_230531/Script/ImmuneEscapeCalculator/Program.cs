using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Diagnostics;

namespace ImmuneEscapeCalculator
{
    public class Antibody
    {
        public string Name;
        public Dictionary<string, double> NeutralDic = new Dictionary<string, double>();//variant, IC50
        public Dictionary<string, double> Escape = new Dictionary<string, double>();//pos mut, escape score
        public double MaxEscape = 0;
        public double AverageEscape = 0;
        public string source;//wt BA.1
        public string group;//12 group
    }
    class Program
    {
        static Dictionary<string, double> Dic_DMS_Exp = new Dictionary<string, double>();//e^(bind+expr)
        static Dictionary<string, int> Dic_Mut_Codon = new Dictionary<string, int>();//突变是否可以通过单个碱基突变到达
        static Dictionary<string, Antibody> Dic_Antibody = new Dictionary<string, Antibody>();//name, antibody
        static List<string> variantList = new List<string>();//测了中和能力的变异株
        static List<string> groupList = new List<string>();//12类抗体
        static string AA20 = "ACDEFGHIKLMNPQRSTVWY";
        static string RBD331_531_AA;
        static string RBD331_531_Nuc;
        List<string> AAMutationList = new List<string>();
        static string Workfold = "M://China220701_230531";
        static string referenceGenome = "+";
        static Dictionary<string, string> mimazi = new Dictionary<string, string>();
        static List<string> GeneList = new List<string>();
        static void ReadinAnnoFile(string fold)
        {
            int i, j, k;
            StreamReader read = new StreamReader(fold + "/mimazi.txt");
            string line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                mimazi.Add(line1[0], line1[1]);
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(fold + "/MN908947.fasta");
            line = read.ReadLine();
            line = read.ReadLine();
            referenceGenome += line;
            read.Close();

            read = new StreamReader(fold + "/GeneCDS.txt");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                GeneList.Add(line);
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static string MutationCombineAnnotation(List<string> background)
        {
            List<string> RBDAAList = new List<string>();
            int i, j, k;
            char[] targetGenome = referenceGenome.ToCharArray();
            if(background[0]!="")
            for (i = 0; i < background.Count(); i++)
            {
                char mutnuc = background[i].Substring(background[i].Length - 1, 1)[0];
                    if ("ATCG".Contains(mutnuc))
                    {
                        int mutpos = Convert.ToInt32(background[i].Substring(1, background[i].Length - 2));
                        targetGenome[mutpos] = (char)mutnuc;
                    }
            }
            for (k = 22517; k <= 23185; k+=3)
            {
                int pos = k;
                string geneName = "TBD";
                int start = 0;
                int end = 0;
                for (i = 0; i < GeneList.Count; i++)
                {
                    string[] line1 = GeneList[i].Split('\t');
                    if (Convert.ToInt32(line1[1]) <= pos && pos <= Convert.ToInt32(line1[2]))
                    {
                        geneName = line1[0];
                        start = Convert.ToInt32(line1[1]);
                        end = Convert.ToInt32(line1[2]);
                        break;
                    }
                }

                if ((pos - start) % 3 == 0)//密码子第一位
                {
                    string refcodon = "" + referenceGenome[pos];
                    refcodon += referenceGenome[pos + 1];
                    refcodon += referenceGenome[pos + 2];
                    string mutcodon = new string(targetGenome, pos, 3);
                    string mut = Convert.ToString((pos - start) / 3 + 1) + mimazi[mutcodon];
                    if (mimazi[refcodon] != mimazi[mutcodon] && !RBDAAList.Contains(mut)) RBDAAList.Add(mut);
                }
                if ((pos - start) % 3 == 1)//密码子第2位
                {
                    string refcodon = "" + referenceGenome[pos - 1];
                    refcodon += referenceGenome[pos];
                    refcodon += referenceGenome[pos + 1];
                    string mutcodon = new string(targetGenome, pos - 1, 3);
                    string mut = Convert.ToString((pos - start) / 3 + 1) + mimazi[mutcodon];
                    if (mimazi[refcodon] != mimazi[mutcodon] && !RBDAAList.Contains(mut)) RBDAAList.Add(mut);
                }
                if ((pos - start) % 3 == 2)//密码子第3位
                {
                    string refcodon = "" + referenceGenome[pos - 2];
                    refcodon += referenceGenome[pos - 1];
                    refcodon += referenceGenome[pos];
                    string mutcodon = new string(targetGenome, pos - 2, 3);
                    string mut = Convert.ToString((pos - start) / 3 + 1) + mimazi[mutcodon];
                    if (mimazi[refcodon] != mimazi[mutcodon] && !RBDAAList.Contains(mut)) RBDAAList.Add(mut);
                }
            }

            string mutline = "NA";
            if(RBDAAList.Count > 0)mutline = RBDAAList[0];
            for (j = 1; j < RBDAAList.Count; j++)
                mutline += "," + RBDAAList[j];
            return mutline;
        }
        static void Readin()//读入数据
        {
            int i, j, k;

            //读入DMS
            StreamReader read = new StreamReader(Workfold + "/Data/final_variant_scores.csv");
            string line;
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split(',');
                string AAmut = line1[2] + line1[3].Substring(0, 1);
                double bindExpr = 0;
                if (line1[6] != "NA") bindExpr += Convert.ToDouble(line1[6]);
                if (line1[13] != "NA") bindExpr += Convert.ToDouble(line1[13]);
                if (line1[0] == "Omicron_BA2")
                    //if(bindExpr >= 0)//【选择突变株】
                    Dic_DMS_Exp.Add(AAmut, bindExpr);//【取这俩加和的e幂】
                //else
                //    Dic_DMS_Exp.Add(AAmut, 1);
                line = read.ReadLine();
            }
            read.Close();


            //读入抗体中和
            //read = new StreamReader(Workfold + "/Data/NeutralWTBA125BF7_Cross.txt");
            read = new StreamReader(Workfold + "/Data/SARS-CoV-2-reinfection-DMS-main/antibody_info.csv");
            line = read.ReadLine();
            string[] linev = line.Split('\t');
            for (i = 1; i <= 15; i++)
                variantList.Add(linev[i]);
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                Antibody newa = new Antibody();
                newa.Name = line1[0];
                double IC50 = 1;
                for (i = 1; i <= 15; i++)
                    if (line1[i] != "NA" && line1[i] != "--")
                    {
                        IC50 = Convert.ToDouble(line1[i]);
                        //if (IC50 > 1) IC50 = 1;
                        //if (IC50 < 0.0005) IC50 = 0.0005;
                        if (IC50 > 10) IC50 = 10;
                        newa.NeutralDic.Add(variantList[i - 1], IC50);
                    }
                    else
                        newa.NeutralDic.Add(variantList[i - 1], 10);
                newa.source = line1[18];
                newa.group = line1[19];
                if (!groupList.Contains(newa.group)) groupList.Add(newa.group);
                if (!Dic_Antibody.ContainsKey(line1[0]))
                    Dic_Antibody.Add(line1[0], newa);

                //存两个键值
                //if(line1[0]!=line1[1])
                //    Dic_Antibody.Add(line1[1], newa);
                line = read.ReadLine();
            }
            read.Close();

            //读入逃逸得分
            read = new StreamReader(Workfold + "/Data/use_res_clean_BA5.txt");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split(',');
                string AAmut = line1[1] + line1[2];
                double Escape = Convert.ToDouble(line1[3]);
                if (Dic_Antibody.ContainsKey(line1[0]))
                {
                    Dic_Antibody[line1[0]].Escape.Add(AAmut, Escape);
                    Dic_Antibody[line1[0]].AverageEscape += Escape;
                    if (Escape > Dic_Antibody[line1[0]].MaxEscape) Dic_Antibody[line1[0]].MaxEscape = Escape;
                }
                line = read.ReadLine();
            }
            read.Close();

            //读入RBD序列，计算密码子
            read = new StreamReader(Workfold + "/Data/RBD_331_531.txt");
            RBD331_531_AA = read.ReadLine();
            RBD331_531_Nuc = read.ReadLine();
            read.Close();
            Dictionary<string, string> CodonMap = new Dictionary<string, string>();
            read = new StreamReader(Workfold + "/Data/mimazi.txt");
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                CodonMap.Add(line1[0], line1[1]);
                line = read.ReadLine();
            }
            read.Close();
            string ATCG = "ATCG";
            for (i = 0; i < 201; i++)
            {
                string refCodon = RBD331_531_Nuc.Substring(i * 3, 3);
                for (j = 0; j < 3; j++)
                {
                    for (k = 0; k < 4; k++)
                    {
                        char[] mutCodon = refCodon.ToCharArray();
                        mutCodon[j] = ATCG[k];
                        string tmps = Convert.ToString(i + 331) + CodonMap[new string(mutCodon)];
                        if (!Dic_Mut_Codon.ContainsKey(tmps))
                            Dic_Mut_Codon.Add(tmps, 1);
                    }
                }
            }

            StreamWriter writeDMS = new StreamWriter(Workfold + "/Data/bindExpr.txt");
            for (i = 331; i <= 531; i++)
            {
                double sum = 0;
                int count = 0;
                for (j = 0; j < AA20.Length; j++)
                {
                    string AAmut = Convert.ToString(i) + Convert.ToString(AA20[j]);
                    double S_dms = 1;
                    if (Dic_DMS_Exp.ContainsKey(AAmut))
                    {
                        S_dms = Dic_DMS_Exp[AAmut];
                        writeDMS.WriteLine(AAmut + "\t" + Dic_DMS_Exp[AAmut]);
                    }
                    if (Dic_Mut_Codon.ContainsKey(AAmut))
                    {
                        sum += S_dms;
                        count++;
                    }

                }
                //writeDMS.WriteLine(Convert.ToString(i) + "\t" + Convert.ToString(sum / count));

            }
            writeDMS.Close();
        }
        static List<double> EscapeScoreCalculator(string Variant, string Antibody)
        {
            int i, j, k;
            List<double> siteScoreList = new List<double>();
            for (i = 331; i <= 531; i++)
            {
                for (j = 0; j < AA20.Length; j++)
                {
                    string AAmut = Convert.ToString(i) + Convert.ToString(AA20[j]);

                    double S_neutral = 0;
                    if (Dic_Antibody[Antibody].NeutralDic[Variant] != 0)
                        S_neutral = -1 * Math.Log10(Dic_Antibody[Antibody].NeutralDic[Variant]);

                    double S_escape = 0;
                    if (Dic_Antibody[Antibody].MaxEscape != 0 && Dic_Antibody[Antibody].Escape.ContainsKey(AAmut))
                        S_escape = Dic_Antibody[Antibody].Escape[AAmut] / Dic_Antibody[Antibody].MaxEscape;

                    double S = S_neutral * S_escape;
                    //double S = S_escape;

                    siteScoreList.Add(S);
                }
            }
            return siteScoreList;
        }
        static void AllAntibodyEscapeScoreCalculator()
        {
            string variantTEsc = "NEW_BF.7";
            StreamWriter write = new StreamWriter(Workfold + "/Data/EscapeScore_PKU_" + variantTEsc + ".single.txt");
            int i, j, k, tmpi, tmpj;
            List<double> TotalEscapeScore = new List<double>();
            List<string> AAmut = new List<string>();
            //先算所有抗体的得分
            for (i = 331; i <= 531; i++)
                for (j = 0; j < AA20.Length; j++)
                {
                    TotalEscapeScore.Add(0);
                    AAmut.Add(Convert.ToString(i) + AA20[j]);
                }
            foreach (string val in Dic_Antibody.Keys)
            {
                if (Dic_Antibody[val].source.Contains("BF.7") || Dic_Antibody[val].source.Contains("reinfection")) //variantTEsc
                {
                    List<double> tmpl = new List<double>(EscapeScoreCalculator("BF7_IC50", val));
                    for (i = 0; i < tmpl.Count; i++)
                        TotalEscapeScore[i] += tmpl[i];
                }
            }
            write.WriteLine("Pos\t" + variantTEsc);
            for (i = 0; i < TotalEscapeScore.Count; i++)
            {
                string output = AAmut[i] + "\t";
                output += Convert.ToString(TotalEscapeScore[i]);
                write.WriteLine(output);
            }
            write.Close();

            //再算每一类的
            groupList.Sort();
            write = new StreamWriter(Workfold + "/Data/EscapeScore_PKU_" + variantTEsc + ".12.txt");
            string outputline = "Mut";
            for (i = 0; i < groupList.Count; i++) outputline += "\t" + groupList[i];
            write.WriteLine(outputline);
            List<string> outgroup = new List<string>();
            for (i = 331; i <= 531; i++)
                for (j = 0; j < AA20.Length; j++)
                {
                    outgroup.Add(Convert.ToString(i) + AA20[j]);
                }
            for (k = 0; k < groupList.Count; k++)
            {
                TotalEscapeScore = new List<double>();
                for (tmpi = 331; tmpi <= 531; tmpi++)
                    for (tmpj = 0; tmpj < AA20.Length; tmpj++)
                        TotalEscapeScore.Add(0);
                foreach (string val in Dic_Antibody.Keys)
                {
                    if (Dic_Antibody[val].group == groupList[k])
                        if (Dic_Antibody[val].source.Contains("BF.7") || Dic_Antibody[val].source.Contains("reinfection"))
                        {
                            List<double> tmpl = new List<double>(EscapeScoreCalculator("BF7_IC50", val));
                            for (tmpi = 0; tmpi < tmpl.Count; tmpi++)
                                if (tmpl[tmpi] > 0)
                                    TotalEscapeScore[tmpi] += tmpl[tmpi];
                        }
                }
                for (i = 0; i < TotalEscapeScore.Count; i++)
                    outgroup[i] += "\t" + Convert.ToString(TotalEscapeScore[i]);
            }
            for (i = 0; i < outgroup.Count; i++)
                write.WriteLine(outgroup[i]);
            write.Close();
        }
        static void EscapeCalculator(string mutdateFile, string variant, string source) //JBloom ImmuneEscapeCalculator
        {
            int i, j, k;
            StreamReader read = new StreamReader(mutdateFile);
            StreamWriter write = new StreamWriter(mutdateFile + ".Pressure");
            //string variant = "BA5_IC50";
            //string source = "BA.5 convalescents";
            string[] filename = mutdateFile.Split('/');
            string[] filename1 = filename[filename.Length - 1].Split('.');
            write.WriteLine("Seq\tDate\tNeutralPressure\tLineage\tRBDCount");
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                string output = line1[0] + "\t" + line1[1] + "\t";
                string RBDAA = MutationCombineAnnotation(line1[2].Split(',').ToList());
                //output += RBDAA + "\t";
                List<string> RBDmuts = RBDAA.Split(',').ToList();
                int A = 0;//TotalAntibodyNumber
                double SumbaM = 0;
                bool RBDChange = false;
                foreach(string key in Dic_Antibody.Keys)
                {
                    if(Dic_Antibody[key].source == source)
                    {
                        A++;
                        double Xar = 1;
                        for (i = 0; i < RBDmuts.Count; i++)
                            if(Dic_Antibody[key].Escape.ContainsKey(RBDmuts[i]))
                                Xar *= (1 - (Dic_Antibody[key].Escape[RBDmuts[i]] / Dic_Antibody[key].MaxEscape));
                            //else
                            //    Xar *= (1 - (Dic_Antibody[key].AverageEscape / Dic_Antibody[key].Escape.Count()) / Dic_Antibody[key].MaxEscape);
                        double Wa = -1 * Math.Log10(Dic_Antibody[key].NeutralDic[variant]);
                        if (Wa < 0) Wa = 0;
                        double baM = Wa * Xar * Xar;
                        SumbaM += baM;
                        if (Xar != 1)
                            RBDChange = true;
                    }
                }
                output += Convert.ToString(SumbaM);
                output += "\t" + filename1[0];
                output += "\t" + Convert.ToString(RBDmuts.Count);
                write.WriteLine(output);
                line = read.ReadLine();
            }
            read.Close();
            write.Close();
            return;
        }
        static void EscapeCalculatorDateMax(string mutdateFile, string variant, string source) //JBloom 逃逸计算器 只记录每个日期的最大值
        {
            int i, j, k;
            StreamReader read = new StreamReader(mutdateFile);
            StreamWriter write = new StreamWriter(mutdateFile + ".MaxPressure");
            StreamWriter writeSeq = new StreamWriter(mutdateFile + ".MaxSequence");
            StreamWriter writeAve = new StreamWriter(mutdateFile + ".AveragePressure");
            Dictionary<string, double> DateMaxEscDic = new Dictionary<string, double>();
            Dictionary<string, string> DateMaxSeqDic = new Dictionary<string, string>();
            Dictionary<string, double> DateAverageEscDic = new Dictionary<string, double>();
            Dictionary<string, double> DateCountEscDic = new Dictionary<string, double>();
            //string variant = "BA5_IC50";
            //string source = "BA.5 convalescents";
            string[] filename = mutdateFile.Split('/');
            string[] filename1 = filename[filename.Length - 1].Split('.');
            write.WriteLine("Seq\tDate\tNeutralPressure\tLineage");
            writeAve.WriteLine("Seq\tDate\tNeutralPressure\tLineage");
            string line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                string output = line1[0] + "\t" + line1[1] + "\t";
                string RBDAA = MutationCombineAnnotation(line1[2].Split(',').ToList());
                //output += RBDAA + "\t";
                List<string> RBDmuts = RBDAA.Split(',').ToList();
                int A = 0;//TotalAntibodyNumber
                double SumbaM = 0;
                bool RBDChange = false;
                foreach (string key in Dic_Antibody.Keys)
                {
                    if (Dic_Antibody[key].source == source)
                    {
                        A++;
                        double Xar = 1;
                        for (i = 0; i < RBDmuts.Count; i++)
                            if (Dic_Antibody[key].Escape.ContainsKey(RBDmuts[i]))
                                Xar *= (1 - (Dic_Antibody[key].Escape[RBDmuts[i]] / Dic_Antibody[key].MaxEscape));
                        //else
                        //    Xar *= (1 - (Dic_Antibody[key].AverageEscape / Dic_Antibody[key].Escape.Count()) / Dic_Antibody[key].MaxEscape);
                        double Wa = -1 * Math.Log10(Dic_Antibody[key].NeutralDic[variant]);
                        if (Wa < 0) Wa = 0;
                        double baM = Wa * Xar * Xar;
                        SumbaM += baM;
                        if (Xar != 1)
                            RBDChange = true;
                    }
                }
                string[] date = line1[1].Split('-');
                if (date.Length > 1)
                {
                    string month = date[0] + "-" + date[1];
                    if (DateMaxEscDic.ContainsKey(month))
                    {
                        if (DateMaxEscDic[month] > SumbaM)
                        {
                            DateMaxEscDic[month] = SumbaM;
                            DateMaxSeqDic[month] = line;
                        }
                        DateAverageEscDic[month] += SumbaM;
                        DateCountEscDic[month] += 1;
                    }
                    else
                    {
                        DateMaxEscDic.Add(month, SumbaM);
                        DateMaxSeqDic.Add(month, line);
                        DateAverageEscDic.Add(month, SumbaM);
                        DateCountEscDic.Add(month, 1);
                    }
                    
                }
                line = read.ReadLine();
            }
            foreach(string key in DateMaxEscDic.Keys)
            {
                string output = key + "\t" + key + "\t";
                output += Convert.ToString(DateMaxEscDic[key]);
                output += "\t" + filename1[0] + "MaxMonth";
                write.WriteLine(output);
                writeSeq.WriteLine(DateMaxSeqDic[key]);

                output = key + "\t" + key + "\t";
                output += Convert.ToString(DateAverageEscDic[key] / DateCountEscDic[key]);
                output += "\t" + filename1[0] + "AverageMonth";
                writeAve.WriteLine(output);
            }

            
            read.Close();
            write.Close();
            writeAve.Close();
            return;
        }
        static void Main(string[] args)
        {
            ReadinAnnoFile(Workfold);
            Readin();//读入数据
            //AllAntibodyEscapeScoreCalculator();
            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/CHN_BF7.MutDate.txt", "BF7_IC50", "BF.7 convalescents");//"BA5_IC50", "BA.5 convalescents"
            EscapeCalculatorDateMax("M://China220701_230531/ChinaVSAbroad/EscVSXBB/CHN_BF7.MutDate.txt", "BF7_IC50", "BF.7 convalescents");
            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BF7.MutDate.txt", "BF7_IC50", "BF.7 convalescents");
            EscapeCalculatorDateMax("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BF7.MutDate.txt", "BF7_IC50", "BF.7 convalescents");
            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BF7_XBB.MutDate.txt", "XBB_IC50", "BF.7 convalescents");
            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BF7_XBB1_5.MutDate.txt", "XBB1_5_IC50", "BF.7 convalescents");
            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BF7_XBB1_5_10.MutDate.txt", "XBB1_5_10_IC50", "BF.7 convalescents");

            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/CHN_BA52.MutDate.txt", "BA5_IC50", "BA.5 convalescents");//"BA5_IC50", "BA.5 convalescents"
            EscapeCalculatorDateMax("M://China220701_230531/ChinaVSAbroad/EscVSXBB/CHN_BA52.MutDate.txt", "BA5_IC50", "BA.5 convalescents");
            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BA52.MutDate.txt", "BA5_IC50", "BA.5 convalescents");
            EscapeCalculatorDateMax("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BA52.MutDate.txt", "BA5_IC50", "BA.5 convalescents");
            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BA5_XBB.MutDate.txt", "XBB_IC50", "BA.5 convalescents");
            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BA5_XBB1_5.MutDate.txt", "XBB1_5_IC50", "BA.5 convalescents");
            EscapeCalculator("M://China220701_230531/ChinaVSAbroad/EscVSXBB/BA5_XBB1_5_10.MutDate.txt", "XBB1_5_10_IC50", "BA.5 convalescents");

            return;
        }
    }
}
