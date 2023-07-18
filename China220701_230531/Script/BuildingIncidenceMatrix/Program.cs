using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace BuildingIncidenceMatrix
{
    public class Lineage
    {
        public string LineageName;
        public Dictionary<string, string> MutIncidenceDic = new Dictionary<string, string>();
        public Dictionary<string, string> RBDAAMutIncidenceDic = new Dictionary<string, string>();
        public List<int> RBDAAMutIncidenceList = new List<int>();
        public double totalRBDMuteventNumber = 0;
    }
    class Program
    {
        static void Main(string[] args)
        {
            string Workfold = "M://China220701_230531";
            List<string> Lineage1 = new List<string>();
            List<string> Lineage2 = new List<string>();
            int i, j, k;
            List<Lineage> LineageList = new List<Lineage>();
            List<string> MutationList = new List<string>();
            List<string> GeneAnnotation = new List<string>();

            List<string> MutationList_NS = new List<string>();
            List<string> GeneAnnotation_NS = new List<string>();

            List<string> MutationList_RBDAA = new List<string>();
            List<string> GeneAnnotation_AARBD = new List<string>();

            string TargetLineage = "BA.5.2";//"BA.5.2" "BF.7"
            
            //Read in sublineage
            List<string> SublineageList = new List<string>();
            StreamReader read = new StreamReader(Workfold + "/ChinaVSAbroad/incidenceCorr/" + TargetLineage + "SublineageList.txt");
            string line = read.ReadLine();
            while(line!=null)
            {
                SublineageList.Add(line);
                line = read.ReadLine();
            }
            read.Close();

            //Read in mutations from the target lineage
            Lineage newl = new Lineage();
            newl.LineageName = TargetLineage;
            read = new StreamReader(Workfold + "/CHNSubtree/MutIncidence/" + TargetLineage + ".MutIncidence");
            line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if (!MutationList.Contains(line1[4]) && line1[6]!="1")//All Mut
                {
                    MutationList.Add(line1[4]);
                    string[] line2 = line1[5].Split('_');
                    if (line2[0].Contains("nsp"))
                        GeneAnnotation.Add("ORF1ab");
                    else if (line2[0] == "S")
                    {
                        int pos = Convert.ToInt32(line2[1].Substring(1, line2[1].Length - 2));
                        if (pos >= 331 && pos <= 531)
                            GeneAnnotation.Add("S_RBD");
                        else
                            if (pos >= 13 && pos <= 305)
                            GeneAnnotation.Add("S_NTD");
                        else
                            GeneAnnotation.Add("S_Other");
                    }
                    else
                        GeneAnnotation.Add(line2[0]);
                }

                string[] linez = line1[5].Split('_');
                if (linez.Length>1 && linez[1][0]!=linez[1][linez[1].Length-1] && !MutationList_NS.Contains(line1[4]) && line1[6] != "1")//NS Mut
                {
                    MutationList_NS.Add(line1[4]);
                    string[] line2 = line1[5].Split('_');
                    if (line2[0].Contains("nsp"))
                        GeneAnnotation_NS.Add("ORF1ab");
                    else if (line2[0] == "S")
                    {
                        int pos = Convert.ToInt32(line2[1].Substring(1, line2[1].Length - 2));
                        if (pos >= 331 && pos <= 531)
                            GeneAnnotation_NS.Add("S_RBD");
                        else
                            if (pos >= 13 && pos <= 305)
                            GeneAnnotation_NS.Add("S_NTD");
                        else
                            GeneAnnotation_NS.Add("S_Other");
                    }
                    else
                        GeneAnnotation_NS.Add(line2[0]);
                }

                newl.MutIncidenceDic.Add(line1[4], line1[6]);
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(Workfold + "/CHNSubtree/MutIncidence/" + TargetLineage + ".MutIncidence.RBDAA");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                if (!MutationList_RBDAA.Contains(line1[0]))//All Mut
                {
                    MutationList_RBDAA.Add(line1[0]);
                }

                newl.RBDAAMutIncidenceDic.Add(line1[0], line1[1]);
                newl.RBDAAMutIncidenceList.Add(Convert.ToInt32(line1[1]));
                newl.totalRBDMuteventNumber += Convert.ToInt32(line1[1]);
                line = read.ReadLine();
            }
            newl.RBDAAMutIncidenceList.Sort();
            LineageList.Add(newl);
            read.Close();

            //Read in mutations from the sub lineage
            for (i=0;i<SublineageList.Count;i++)
            {
                read = new StreamReader(Workfold + "/CHNSubtree/MutIncidence/" + SublineageList[i] + ".MutIncidence");
                Lineage newL = new Lineage();
                newL.LineageName = SublineageList[i];
                line = read.ReadLine();
                line = read.ReadLine();
                while (line != null)
                {
                    string[] line1 = line.Split('\t');
                    string[] line2 = line1[5].Split('_');
                    if (!MutationList.Contains(line1[4]))
                    {
                        MutationList.Add(line1[4]);
                        if (line2[0].Contains("nsp"))
                            GeneAnnotation.Add("ORF1ab");
                        else if (line2[0] == "S")
                        {
                            int pos = Convert.ToInt32(line2[1].Substring(1, line2[1].Length - 2));
                            if (pos >= 331 && pos <= 531)
                                GeneAnnotation.Add("S_RBD");
                            else
                                if (pos >= 13 && pos <= 305)
                                GeneAnnotation.Add("S_NTD");
                            else
                                GeneAnnotation.Add("S_Other");
                        }
                        else
                            GeneAnnotation.Add(line2[0]);
                    }

                    if (line2.Length > 1 && line2[1][0] != line2[1][line2[1].Length - 1] && !MutationList_NS.Contains(line1[4]))//NS Mut
                    {
                        MutationList_NS.Add(line1[4]);
                        if (line2[0].Contains("nsp"))
                            GeneAnnotation_NS.Add("ORF1ab"); 
                        else if (line2[0] == "S")
                        {
                            int pos = Convert.ToInt32(line2[1].Substring(1, line2[1].Length - 2));
                            if (pos >= 331 && pos <= 531)
                                GeneAnnotation_NS.Add("S_RBD");
                            else
                                if (pos >= 13 && pos <= 305)
                                GeneAnnotation_NS.Add("S_NTD");
                            else
                                GeneAnnotation_NS.Add("S_Other");
                        }
                        else
                            GeneAnnotation_NS.Add(line2[0]);
                    }

                    if (newL.MutIncidenceDic.ContainsKey(line1[4]))
                    {
                        newL.MutIncidenceDic[line1[4]] = Convert.ToString(Convert.ToInt32(newL.MutIncidenceDic[line1[4]]) + Convert.ToInt32(line1[6]));
                    }
                    else
                        newL.MutIncidenceDic.Add(line1[4], line1[6]);
                    line = read.ReadLine();
                }

                read = new StreamReader(Workfold + "/CHNSubtree/MutIncidence/" + SublineageList[i] + ".MutIncidence.RBDAA");
                line = read.ReadLine();
                line = read.ReadLine();
                while (line != null)
                {
                    string[] line1 = line.Split('\t');
                    if (!MutationList_RBDAA.Contains(line1[0]))//All Mut
                    {
                        MutationList_RBDAA.Add(line1[0]);
                    }

                    newL.RBDAAMutIncidenceDic.Add(line1[0], line1[1]);
                    newL.RBDAAMutIncidenceList.Add(Convert.ToInt32(line1[1]));
                    newL.totalRBDMuteventNumber += Convert.ToInt32(line1[1]);
                    line = read.ReadLine();
                }
                newL.RBDAAMutIncidenceList.Sort();
                LineageList.Add(newL);
                read.Close();
            }

            //output
            StreamWriter write = new StreamWriter(Workfold + "/ChinaVSAbroad/incidenceCorr/" + TargetLineage + "SublineageMuts.txt");
            StreamWriter writeNS = new StreamWriter(Workfold + "/ChinaVSAbroad/incidenceCorr/" + TargetLineage + "SublineageNSMuts.txt");
            StreamWriter writeRBDAA = new StreamWriter(Workfold + "/ChinaVSAbroad/incidenceCorr/" + TargetLineage + "SublineageRBDAAMuts.txt");
            StreamWriter wrtieRBDAATopN = new StreamWriter(Workfold + "/ChinaVSAbroad/incidenceCorr/" + TargetLineage + "SublineageRBDAATopMuts.txt");
            string title = "Gene\tMutation";
            for (j = 0; j < LineageList.Count; j++)
                title += "\t" + LineageList[j].LineageName;
            write.WriteLine(title);
            writeNS.WriteLine(title);
            title = "Mutation";
            for (j = 0; j < LineageList.Count; j++)
                title += "\t" + LineageList[j].LineageName;
            writeRBDAA.WriteLine(title);
            wrtieRBDAATopN.WriteLine(title);

            for (i = 0; i < MutationList.Count; i++)
            {
                string output = GeneAnnotation[i] + "\t" + MutationList[i];
                for (j = 0; j < LineageList.Count; j++)
                {
                    if (LineageList[j].MutIncidenceDic.ContainsKey(MutationList[i]))
                        output += "\t" + LineageList[j].MutIncidenceDic[MutationList[i]];
                    else
                        output += "\t0";//这里看缺失值数据怎么处理
                }
                write.WriteLine(output);
            }
            for (i = 0; i < MutationList_NS.Count; i++)
            {
                string output = GeneAnnotation_NS[i] + "\t" + MutationList_NS[i];
                for (j = 0; j < LineageList.Count; j++)
                {
                    if (LineageList[j].MutIncidenceDic.ContainsKey(MutationList_NS[i]))
                        output += "\t" + LineageList[j].MutIncidenceDic[MutationList_NS[i]];
                    else
                        output += "\t0";//这里看缺失值数据怎么处理
                }
                writeNS.WriteLine(output);
            }

            List<string> RBDTopN = new List<string>();
            for (i = 0; i < MutationList_RBDAA.Count; i++)
            {
                string output = MutationList_RBDAA[i];
                for (j = 0; j < LineageList.Count; j++)
                {
                    if (LineageList[j].RBDAAMutIncidenceDic.ContainsKey(MutationList_RBDAA[i]))
                    {
                        output += "\t" + Convert.ToString(Convert.ToDouble(LineageList[j].RBDAAMutIncidenceDic[MutationList_RBDAA[i]]) / LineageList[j].totalRBDMuteventNumber);
                        if (Convert.ToDouble(LineageList[j].RBDAAMutIncidenceDic[MutationList_RBDAA[i]]) >= LineageList[j].RBDAAMutIncidenceList[LineageList[j].RBDAAMutIncidenceList.Count - 5])
                            if (!RBDTopN.Contains(MutationList_RBDAA[i]))
                                RBDTopN.Add(MutationList_RBDAA[i]);
                    }
                    else
                        output += "\t0";//这里看缺失值数据怎么处理
                    
                }
                writeRBDAA.WriteLine(output);
            }

            for (i = 0; i < RBDTopN.Count; i++)
            {
                string output = RBDTopN[i];
                for (j = 0; j < LineageList.Count; j++)
                {
                    if (LineageList[j].RBDAAMutIncidenceDic.ContainsKey(RBDTopN[i]))
                    {
                        output += "\t" + Convert.ToString(Convert.ToDouble(LineageList[j].RBDAAMutIncidenceDic[RBDTopN[i]]) / LineageList[j].totalRBDMuteventNumber);
                    }
                    else
                        output += "\t0";//这里看缺失值数据怎么处理
                }
                wrtieRBDAATopN.WriteLine(output);
            }
            write.Close();
            writeNS.Close();
            writeRBDAA.Close();
            wrtieRBDAATopN.Close();
        }
    }
}
