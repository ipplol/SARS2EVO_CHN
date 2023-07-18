using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace FindCHNSubtreeFromMAT
{
    public class Node
    {
        public string name;
        public int FatherInD;
        public bool IsLeaf = true;
        public List<string> Nuc_mutation = new List<string>();
        public List<string> Total_Nuc_mutation = new List<string>();
        public string CollectionDate = null;
        public string Location = "TBD";
        public string Lineage;
        public string NextstrainClade;
        public string accessionID = "NA";
        public List<int> ChildernIndList = new List<int>();//子节点编号
        public int SameSequenceNumber = 0;//和这个节点一样的叶子有多少
        
        public bool NewCHN = false;//这个节点是新增的中国序列吗
        public int Num_AllDescendants = 0;//该节点共有多少子代叶子节点
        public int Num_CHNDescendants = 0;//该节点共有多少中国子代叶子节点
        public int Num_CHNNewDescendants = 0;//该节点共有多少新增的中国子代叶子节点
    }
    public class Metadata
    {
        public string name;
        public string CollectionDate;
        public string Location;
        public string Lineage;
        public string NextstrainClade = "Other";
        public string Accession = "TBD";
    }
    class Program
    {
        static List<string> JsonLine = new List<string>(234217729);
        static List<Node> NodesList = new List<Node>();
        static Dictionary<string, Metadata> MetadataDic = new Dictionary<string, Metadata>();
        static Dictionary<string, Metadata> CHNSeqMetaDic = new Dictionary<string, Metadata>();
        static string JFile = "";
        static string Workfold = "";
        static string referenceGenome = "+";
        static Dictionary<string, string> mimazi = new Dictionary<string, string>();
        static List<string> GeneList = new List<string>();

        static void ReadinAnnoFile(string fold)
        {
            int i, j, k;
            StreamReader read = new StreamReader(fold + "/data/mimazi.txt");
            string line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                mimazi.Add(line1[0], line1[1]);
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(fold + "/data/MN908947.fasta");
            line = read.ReadLine();
            line = read.ReadLine();
            referenceGenome += line;
            read.Close();

            read = new StreamReader(fold + "/data/GeneCDS.txt");
            line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                GeneList.Add(line);
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static string ColorPicker(int nPicks)//选择n个差异最大的16进制颜色
        {
            // Define the list of hexadecimal RGB codes
            //List<string> hexColors = new List<string> { "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF", "#800000", "#008000", "#000080" };
            List<string> hexColors = new List<string> ("#990033,#CC6699,#FF6699,#FF3366,#993366,#CC0066,#CC0033,#FF0066,#FF3399,#FF9999,#FF99CC,#FF0099,#CC3366,#FF66CC,#FF33CC,#FFCCFF,#FF0033,#FF00CC,#CC3399,#FF99FF,#FF66FF,#CC33CC,#CC00FF,#FF33FF,#CC99FF,#9900CC,#FF00FF,#CC66FF,#CC33FF,#CC99CC,#990066,#993399,#CC66CC,#CC00CC,#663366,#CC0099,#990099,#660099,#666FFF,#0000CC,#9933CC,#666699,#660066,#333366,#0066CC,#99CCFF,#9933FF,#330099,#6699FF,#9966CC,#3300CC,#003366,#330033,#663399,#3333FF,#006699,#6633CC,#3333CC,#3399CC,#6600CC,#0066FF,#0033FF,#66CCFF,#330066,#3366FF,#3399FF,#6600FF,#3366CC,#6699CC,#0099FF,#CCCCFF,#000033,#33CCFF,#9999FF,#0000FF,#00CCFF,#9999CC,#0033CC,#3300FF,#333399,#000099,#000066,#6633FF,#003399,#6666CC,#0099CC,#9900FF,#9966FF,#FFFFCC,#FFCC00,#CC9909,#663300,#FF6600,#663333,#CC6666,#FF6666,#FFCC66,#FF9900,#FF9966,#CC3300,#996666,#FFCCCC,#660000,#FF3300,#CC6600,#FF6633,#996633,#CC9999,#FF3333,#990000,#CC9966,#FFFF33,#FF9933,#330000,#993333,#CC3333,#CC0000,#FFCC99,#FFFF00,#996600,#993300,#FF0000,#CC6633,#CC9933,#FFCC33,#FFFF99,#99FFFF,#33CCCC,#00CC99,#99FF99,#009966,#33FF33,#33FF00,#99CC33,#66CCCC,#66FFCC,#66FF66,#009933,#00CC33,#66FF00,#336600,#333000,#99FFCC,#339933,#33FF66,#33CC33,#99FF00,#669900,#666600,#00FFFF,#99CC99,#00FF66,#66FF33,#66CC00,#99CC00,#999933,#00CCCC,#006666,#CCFFCC,#00FF00,#00CC00,#CCFF66,#CCCC66,#009999,#003333,#006633,#66CC33,#33CC00,#CCFF33,#666633,#669999,#00FFCC,#336633,#33CC66,#339900,#CCFF00,#999966,#99CCCC,#33FFCC,#669966,#00CC66,#99FF33,#999900,#CCCC99,#CCFFFF,#33CC99,#66CC66,#66CC99,#00FF33,#009900,#CCCC00,#CCC333,#336666,#006600,#003300,#669933,#339966,#339999,#669900,#99CC66,#99FF66,#00FF99,#33FF99,#66FF99,#CCFF99,#33FFFF,#66FFFF".Split(','));
            // Convert each hexadecimal RGB code to its corresponding RGB value
            List<(int R, int G, int B)> rgbColors = new List<(int R, int G, int B)>();
            foreach (string hexColor in hexColors)
            {
                int r = Convert.ToInt32(hexColor.Substring(1, 2), 16);
                int g = Convert.ToInt32(hexColor.Substring(3, 2), 16);
                int b = Convert.ToInt32(hexColor.Substring(5, 2), 16);
                rgbColors.Add((r, g, b));
            }

            // Calculate the pairwise distances between all colors
            int nColors = hexColors.Count;
            double[,] distances = new double[nColors, nColors];
            for (int i = 0; i < nColors; i++)
            {
                for (int j = i + 1; j < nColors; j++)
                {
                    distances[i, j] = Math.Sqrt(Math.Pow(rgbColors[i].R - rgbColors[j].R, 2) + Math.Pow(rgbColors[i].G - rgbColors[j].G, 2) + Math.Pow(rgbColors[i].B - rgbColors[j].B, 2));
                    distances[j, i] = distances[i, j];
                }
            }

            // Pick N colors that are the most different from each other
            List<string> pickedColors = new List<string>();
            pickedColors.Add(hexColors[nPicks]); // start with the first color
            for (int i = 1; i < nPicks; i++)
            {
                double maxMinDistance = double.MinValue;
                string maxMinDistanceColor = "";
                foreach (string hexColor in hexColors)
                {
                    if (pickedColors.Contains(hexColor))
                    {
                        continue;
                    }
                    double minDistance = double.MaxValue;
                    foreach (string pickedColor in pickedColors)
                    {
                        int index1 = hexColors.IndexOf(hexColor);
                        int index2 = hexColors.IndexOf(pickedColor);
                        double distance = distances[index1, index2];
                        if (distance < minDistance)
                        {
                            minDistance = distance;
                        }
                    }
                    if (minDistance > maxMinDistance)
                    {
                        maxMinDistance = minDistance;
                        maxMinDistanceColor = hexColor;
                    }
                }
                pickedColors.Add(maxMinDistanceColor);
            }

            return string.Join(", ", pickedColors).Replace(" ","");
        }
        static bool DateCompare(string A, string B)//A>B true else false
        {
            int i;
            for (i = 0; i < A.Length && i < B.Length; i++)
            {
                if (A[i] > B[i])
                    return true;
            }
            return false;
        }
        static List<string> FindEarliestDate(List<string> DateList, List<string> LocationList)//从list中找出最早采样时间
        {
            string tmps;
            int i, j;
            for (i = 0; i < DateList.Count; i++)
            {
                for (j = i + 1; j < DateList.Count; j++)
                {
                    if (DateCompare(DateList[i], DateList[j]))
                    {
                        tmps = DateList[i]; DateList[i] = DateList[j]; DateList[j] = tmps;
                        tmps = LocationList[i]; LocationList[i] = LocationList[j]; LocationList[j] = tmps;
                    }
                }
            }
            List<string> RT = new List<string>();
            RT.Add(DateList[0]);
            RT.Add(LocationList[0]);
            return RT;
        }
        static void AssignDateToNodes()//给树中间的节点分配采样时间
        {
            /*
             * 树上存在完全一致的序列，这些序列难以判断发生突变，因为他们到前一个节点的距离都是0
             * 前一个节点，也是树中间节点的突变，可以代表这些序列的突变
             * 用这些序列的最早采样时间，作为前一个节点的采样时间
             */
            int i, j, k;
            for (i = 0; i < NodesList.Count; i++)
            {
                if (!NodesList[i].IsLeaf)
                {
                    List<string> DateList = new List<string>();
                    List<string> LocationList = new List<string>();
                    for (j = 0; j < NodesList[i].ChildernIndList.Count; j++)
                    {
                        //该中间节点有叶子节点的子代且子代没有额外突变
                        if (NodesList[NodesList[i].ChildernIndList[j]].IsLeaf && NodesList[NodesList[i].ChildernIndList[j]].Nuc_mutation.Count == 0)
                            if (NodesList[NodesList[i].ChildernIndList[j]].CollectionDate != null)
                            {
                                DateList.Add(NodesList[NodesList[i].ChildernIndList[j]].CollectionDate);
                                LocationList.Add(NodesList[NodesList[i].ChildernIndList[j]].Location);
                            }
                        NodesList[i].SameSequenceNumber = DateList.Count;
                    }
                    if (DateList.Count > 0)
                    {
                        List<string> DateLocation = FindEarliestDate(DateList, LocationList);
                        NodesList[i].CollectionDate = DateLocation[0];
                        NodesList[i].Location = DateLocation[1];
                        NodesList[i].Lineage = NodesList[NodesList[i].ChildernIndList[0]].Lineage;
                        NodesList[i].NextstrainClade = NodesList[NodesList[i].ChildernIndList[0]].NextstrainClade;
                    }

                }
            }

            return;
        }
        static void LeafFindMetadata()//给叶子节点找metadata
        {
            int i, j, k;
            for (i = 0; i < NodesList.Count; i++)
            {
                if (NodesList[i].IsLeaf)
                {
                    string[] name1 = NodesList[i].name.Split('|');
                    if (CHNSeqMetaDic.ContainsKey(name1[0]))//是不是新加入的中国序列
                    {
                        NodesList[i].accessionID = CHNSeqMetaDic[name1[0]].Accession;
                        NodesList[i].CollectionDate = CHNSeqMetaDic[name1[0]].CollectionDate;
                        if (NodesList[i].CollectionDate.Length == 4)
                            NodesList[i].CollectionDate += "-99-99";
                        if (NodesList[i].CollectionDate.Length == 7)
                            NodesList[i].CollectionDate += "-99";
                        NodesList[i].Lineage = CHNSeqMetaDic[name1[0]].Lineage;
                        NodesList[i].NextstrainClade = CHNSeqMetaDic[name1[0]].NextstrainClade;
                        NodesList[i].Location = CHNSeqMetaDic[name1[0]].Location;
                        NodesList[i].NewCHN = true;
                    }
                    else
                    if (MetadataDic.ContainsKey(NodesList[i].name))
                    {
                        NodesList[i].accessionID = NodesList[i].name;
                        NodesList[i].CollectionDate = MetadataDic[NodesList[i].name].CollectionDate;
                        if (NodesList[i].CollectionDate.Length == 4)
                            NodesList[i].CollectionDate += "-99-99";
                        if (NodesList[i].CollectionDate.Length == 7)
                            NodesList[i].CollectionDate += "-99";
                        NodesList[i].Lineage = MetadataDic[NodesList[i].name].Lineage;
                        NodesList[i].NextstrainClade = MetadataDic[NodesList[i].name].NextstrainClade;
                        NodesList[i].Location = MetadataDic[NodesList[i].name].Location;
                    }
                    else
                        if (name1.Count() > 1 && MetadataDic.ContainsKey(name1[1]))//搜索accession id
                    {
                        NodesList[i].accessionID = name1[1];
                        NodesList[i].CollectionDate = MetadataDic[name1[1]].CollectionDate;
                        if (NodesList[i].CollectionDate.Length == 4)
                            NodesList[i].CollectionDate += "-99-99";
                        if (NodesList[i].CollectionDate.Length == 7)
                            NodesList[i].CollectionDate += "-99";
                        NodesList[i].Lineage = MetadataDic[name1[1]].Lineage;
                        NodesList[i].NextstrainClade = MetadataDic[name1[1]].NextstrainClade;
                        NodesList[i].Location = MetadataDic[name1[1]].Location;
                    }
                    else
                        if (MetadataDic.ContainsKey(name1[0]))//找不到看看文件名有没有
                    {
                        NodesList[i].accessionID = "Seq_Name";
                        NodesList[i].CollectionDate = MetadataDic[name1[0]].CollectionDate;
                        if (NodesList[i].CollectionDate.Length == 4)
                            NodesList[i].CollectionDate += "-99-99";
                        if (NodesList[i].CollectionDate.Length == 7)
                            NodesList[i].CollectionDate += "-99";
                        NodesList[i].Lineage = MetadataDic[name1[0]].Lineage;
                        NodesList[i].NextstrainClade = MetadataDic[name1[0]].NextstrainClade;
                        NodesList[i].Location = MetadataDic[name1[0]].Location;
                    }
                    else//确实找不到，看看文件名里写没写
                    {
                        string[] name11 = name1[0].Split('/');
                        if (name11.Length == 3)
                        {
                            NodesList[i].CollectionDate = name1[name1.Length - 1];
                            if (NodesList[i].CollectionDate.Length == 4)
                                NodesList[i].CollectionDate += "-99-99";
                            if (NodesList[i].CollectionDate.Length == 7)
                                NodesList[i].CollectionDate += "-99";
                            NodesList[i].Location = name11[0];
                        }
                    }
                }
            }
            return;
        }
        static void TreeBuild(string jsonLine, int fatherInd)
        {
            Node newNode = new Node();
            int i, j, k;
            int kuohao = 0;
            //Console.WriteLine(jsonLine);

            //---------分割序列---------
            List<int> SplitPos = new List<int>();//直系逗号的位置
            List<int> ChildSplitPos = new List<int>();//孩子分割位置
            SplitPos.Add(0);
            for (i = 0; i < jsonLine.Length; i++)
            {
                if (kuohao == 0 && jsonLine[i] == ',')
                    SplitPos.Add(i + 1);
                if (jsonLine[i] == '{' || jsonLine[i] == '[') kuohao++;
                if (jsonLine[i] == '}' || jsonLine[i] == ']') kuohao--;
            }
            for (i = 0; i < SplitPos.Count; i++)//切割序列，找突变和子代
            {
                j = SplitPos[i] + 1;
                while (jsonLine[j] != '\"') j++;
                string tag = jsonLine.Substring(SplitPos[i] + 1, j - SplitPos[i] - 1);
                if (tag == "branch_attrs")
                {
                    string lineba1 = jsonLine.Substring(SplitPos[i] + 1, SplitPos[i + 1] - SplitPos[i] - 2);
                    string lineba2 = lineba1.Replace('}', '{');
                    string[] lineba3 = lineba2.Split('{');
                    if (lineba3[4] != "\"nuc\":[]")
                    {
                        string lineba4 = lineba3[4].Substring(7, lineba3[4].Length - 8);
                        string[] lineba5 = lineba4.Split(',');
                        for (k = 0; k < lineba5.Length; k++)
                        {
                            newNode.Nuc_mutation.Add(lineba5[k].Substring(1, lineba5[k].Length - 2));
                        }
                    }
                }
                if (tag == "children")
                {
                    ChildSplitPos.Add(SplitPos[i] + 11);
                    kuohao = 0;
                    for (k = SplitPos[i] + 12; k < SplitPos[i + 1]; k++)
                    {
                        if (kuohao == 0 && jsonLine[k] == ',')
                            ChildSplitPos.Add(k);
                        if (jsonLine[k] == '{' || jsonLine[k] == '[') kuohao++;
                        if (jsonLine[k] == '}' || jsonLine[k] == ']') kuohao--;
                    }
                    ChildSplitPos.Add(SplitPos[i + 1] - 2);
                }
                if (tag == "node_attrs")
                {

                }
                if (tag == "name")
                {
                    string lineba1 = jsonLine.Substring(SplitPos[i] + 1, SplitPos[i + 1] - SplitPos[i] - 2);
                    newNode.name = lineba1.Substring(7, lineba1.Length - 8);
                }
            }

            //---------从父节点继承突变----------
            newNode.Total_Nuc_mutation = new List<string>(newNode.Nuc_mutation);
            if (fatherInd != -1)
            {
                for (i = 0; i < NodesList[fatherInd].Total_Nuc_mutation.Count; i++)
                    newNode.Total_Nuc_mutation.Add(NodesList[fatherInd].Total_Nuc_mutation[i]);
            }

            //--------递归-----------
            if (ChildSplitPos.Count <= 2)//没有子代
            {
                newNode.IsLeaf = true;
                NodesList[fatherInd].ChildernIndList.Add(NodesList.Count);
                NodesList.Add(newNode);
            }
            else//递归子代
            {
                newNode.IsLeaf = false;
                if (fatherInd != -1)
                    NodesList[fatherInd].ChildernIndList.Add(NodesList.Count);
                NodesList.Add(newNode);
                for (k = 0; k < ChildSplitPos.Count - 1; k++)
                {
                    //Console.WriteLine(jsonLine.Substring(ChildSplitPos[k] + 1, ChildSplitPos[k + 1] - ChildSplitPos[k] - 1));
                    TreeBuild(jsonLine.Substring(ChildSplitPos[k] + 2, ChildSplitPos[k + 1] - ChildSplitPos[k] - 3), NodesList.Count - 1);
                }
            }
            return;
        }
        static int TreeBuild2(int lineId, int fatherInd)
        {
            int i, j, k;
            int thisNode = -1;
            int kuohao = 0;
            Node newnode = new Node();
            newnode.FatherInD = fatherInd;
            //继承父节点突变
            if (fatherInd != -1)
                for (k = 0; k < NodesList[fatherInd].Total_Nuc_mutation.Count; k++)
                    newnode.Total_Nuc_mutation.Add(NodesList[fatherInd].Total_Nuc_mutation[k]);
            while (lineId < JsonLine.Count())
            {
                if (JsonLine[lineId].IndexOf('{') != -1) kuohao++;
                if (JsonLine[lineId].IndexOf('[') != -1) kuohao++;
                if (JsonLine[lineId].IndexOf('}') != -1) kuohao--;
                if (JsonLine[lineId].IndexOf(']') != -1) kuohao--;
                if (kuohao == 0)
                    break;
                //------------节点名称-----------
                if (JsonLine[lineId].Contains("nuc mutations"))
                {
                    string[] label = JsonLine[lineId].Split(':');
                    newnode.name = label[3].Substring(1, label[3].Length - 3);
                }
                //------------节点突变-----------
                if (JsonLine[lineId] == "\"nuc\":[")
                {
                    lineId++; kuohao--;
                    string[] mut1 = JsonLine[lineId].Split('\"');
                    for (i = 0; i < mut1.Count(); i++)
                        if (mut1[i].Length > 2)
                        {
                            newnode.Nuc_mutation.Add(mut1[i]);
                            newnode.Total_Nuc_mutation.Add(mut1[i]);
                        }
                }
                //------------递归子代-----------
                if (JsonLine[lineId] == "\"children\":[")
                {
                    newnode.IsLeaf = false;
                    NodesList.Add(newnode);
                    thisNode = NodesList.Count() - 1;
                    lineId++;
                    while (JsonLine[lineId] != "]")
                        lineId = TreeBuild2(lineId, thisNode);
                    if (JsonLine[lineId] == "]") kuohao--;
                }
                lineId++;
            }
            //--------------和父节点建立联系--------------
            if (newnode.IsLeaf)
            {
                NodesList.Add(newnode);
                thisNode = NodesList.Count() - 1;
            }
            if (fatherInd != -1)
            {
                //for (k = 0; k < NodesList[fatherInd].Total_Nuc_mutation.Count; k++)
                //    NodesList[thisNode].Total_Nuc_mutation.Add(NodesList[fatherInd].Total_Nuc_mutation[k]);
                NodesList[fatherInd].ChildernIndList.Add(thisNode);
            }

            return lineId + 1;
        }
        static void MutationEventCal(string jsfile)
        {
            //-----------------------读取json转换输出----------------------
            int i, j, k;
            FileStream fs = new FileStream(jsfile, FileMode.Open);
            TextReader read = new StreamReader(fs);
            StreamWriter writel = new StreamWriter(jsfile + ".line");
            var clen = 1024 * 1024;
            var buffer = new Char[clen];
            var count = read.Read(buffer, 0, clen);
            while (count > 0)
            {
                var str = new string(buffer, 0, count);
                StringBuilder strBuilder = new StringBuilder();
                strBuilder.Append(str);
                strBuilder.Replace("{", "{\n");
                strBuilder.Replace("[", "[\n");
                strBuilder.Replace("}", "}\n");
                strBuilder.Replace("]", "]\n");
                strBuilder.Replace(",", "");
                writel.Write(strBuilder.ToString());
                count = read.Read(buffer, 0, clen);
            }
            read.Close();
            writel.Close();
            read.Dispose();
            Console.WriteLine("Json Reform Done.");

            //------------------------读取新格式json----------------
            StreamReader readline = new StreamReader(jsfile + ".line");
            string line = readline.ReadLine();
            while (line != null)
            {
                JsonLine.Add(line);
                line = readline.ReadLine();
            }
            Console.WriteLine("Read Json Line Done.");
            readline.Close();
            //------------------------构建树结构---------------------
            k = 0;
            while (JsonLine[k] != "\"tree\":{") k++;
            k++;
            TreeBuild2(k, -1);
            Console.WriteLine("Tree Build Done.");

            //-------------------------------------------------


            //---------------------------------------------------------
            LeafFindMetadata();//给叶子节点找metadata
            Console.WriteLine("Finding Metadata.");

            AssignDateToNodes();//给树中间的节点分配采样时间
            Console.WriteLine("Assign Sampling Date.");

            
            //-------------------------输出结果 突变事件------------------------
            Console.WriteLine("Writing Result ......");
            StreamWriter write = new StreamWriter(jsfile + ".mutevent");
            //StreamWriter writeX = new StreamWriter(jsfile + ".XBB");
            write.Write("Sample\tAccessionID\tOriginalMut\tNewMut\tLineage\tCollectionDate\tLocation\tSameSeqNumber\tNextstrainClade\n");
            for (i = 0; i < NodesList.Count; i++)
            {
                //if (NodesList[i].Lineage==null || NodesList[i].Lineage.Contains("XBB"))
                //    writeX.WriteLine(NodesList[i].name + "\t" + Convert.ToString(NodesList[i].Nuc_mutation.Count));
                if ((NodesList[i].IsLeaf || NodesList[i].SameSequenceNumber > 0) && NodesList[i].Nuc_mutation.Count == 1)
                {
                    string output = NodesList[i].name;
                    output += "\t" + NodesList[i].accessionID;
                    string mutation = "";
                    for (j = 0; j < NodesList[i].Total_Nuc_mutation.Count; j++)
                        mutation += NodesList[i].Total_Nuc_mutation[j] + ",";
                    if (mutation != "")
                        mutation = mutation.Substring(0, mutation.Length - 1);
                    output += "\t" + mutation;
                    output += "\t" + NodesList[i].Nuc_mutation[0];
                    output += "\t" + NodesList[i].Lineage;
                    output += "\t" + NodesList[i].CollectionDate;
                    output += "\t" + NodesList[i].Location;
                    output += "\t" + Convert.ToString(NodesList[i].SameSequenceNumber);
                    output += "\t" + NodesList[i].NextstrainClade;
                    write.WriteLine(output);
                }
                else
                    if ((NodesList[i].IsLeaf || NodesList[i].SameSequenceNumber > 0) && NodesList[i].Nuc_mutation.Count == 2)
                {
                    string output = NodesList[i].name;
                    output += "\t" + NodesList[i].accessionID;
                    string mutation = "";
                    for (j = 0; j < NodesList[i].Total_Nuc_mutation.Count; j++)
                        mutation += NodesList[i].Total_Nuc_mutation[j] + ",";
                    if (mutation != "")
                        mutation = mutation.Substring(0, mutation.Length - 1);
                    output += "\t" + mutation;
                    string output1 = output, output2 = output;
                    output1 += "\t" + NodesList[i].Nuc_mutation[0];
                    output1 += "\t" + NodesList[i].Lineage;
                    output1 += "\t" + NodesList[i].CollectionDate;
                    output1 += "\t" + NodesList[i].Location;
                    output1 += "\t" + Convert.ToString(NodesList[i].SameSequenceNumber);
                    output1 += "\t" + NodesList[i].NextstrainClade;
                    write.WriteLine(output1);
                    output2 += "\t" + NodesList[i].Nuc_mutation[1];
                    output2 += "\t" + NodesList[i].Lineage;
                    output2 += "\t" + NodesList[i].CollectionDate;
                    output2 += "\t" + NodesList[i].Location;
                    output2 += "\t" + Convert.ToString(NodesList[i].SameSequenceNumber);
                    output2 += "\t" + NodesList[i].NextstrainClade;
                    write.WriteLine(output2);
                }
            }
            write.Close();
            //writeX.Close();
            return;
        }
        static void ReadMetadata()
        {
            string file;
            StreamReader read;
            string line;
            int i, j, k;
            read = new StreamReader("M://China220701_230531/Data/public-latestCHNadded.metadata.tsv");
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                Metadata newm = new Metadata();
                newm.name = line1[0];
                newm.Lineage = line1[8];
                newm.NextstrainClade = line1[7];
                newm.CollectionDate = line1[2];
                newm.Location = line1[3];
                if (!MetadataDic.ContainsKey(line1[0]))
                    MetadataDic.Add(line1[0], newm);
                if (!MetadataDic.ContainsKey(line1[1]))
                    MetadataDic.Add(line1[1], newm);
                line = read.ReadLine();
            }
            read.Close();
            
            read = new StreamReader("M://China220701_230531/Data/public_allCHN.metadata.tsv");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line!=null)
            {
                string[] line1 = line.Split('\t');
                Metadata newm = new Metadata();
                newm.name = line1[0];
                newm.Lineage = line1[8];
                newm.NextstrainClade = line1[7];
                newm.CollectionDate = line1[2];
                newm.Location = line1[3];
                newm.Accession = line1[1];
                if (!CHNSeqMetaDic.ContainsKey(line1[0]))
                    CHNSeqMetaDic.Add(line1[0], newm);
                line = read.ReadLine();
            }
            read.Close();

            /*read = new StreamReader("M://metadata.tsv");
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                Metadata newm = new Metadata();
                newm.name = line1[0];
                newm.Lineage = line1[4];
                newm.CollectionDate = line1[10];
                string[] line2 = line1[11].Split('/');
                k = line2[0].Length - 1;
                while (line2[0][k] == ' ') k--;
                newm.Location = line2[0].Substring(0, k + 1);
                if (line1[1] != "null" && line1[1] != "" && !MetadataDic.ContainsKey(line1[1]))
                    MetadataDic.Add(line1[1], newm);
                if (line1[3] != "null" && line1[3] != "" && line1[3] != " " && !MetadataDic.ContainsKey(line1[3]))
                    MetadataDic.Add(line1[3], newm);
                if (!MetadataDic.ContainsKey(line1[0]))
                    MetadataDic.Add(line1[0], newm);
                line = read.ReadLine();
            }
            read.Close();*/
            Console.WriteLine("Read Metadata Done");
        }
        static void WriteNodeInfo(string jsfile)
        {
            int i, j, k;
            StreamWriter write = new StreamWriter(jsfile + ".nodeinfo");
            write.WriteLine("NodeName\tEarliestDate\tEarliestLocation\tLineage\tGenomeMut");
            for (i = 0; i < NodesList.Count; i++)
            {
                string output = NodesList[i].name;
                output += "\t" + NodesList[i].CollectionDate;
                output += "\t" + NodesList[i].Location;
                output += "\t" + NodesList[i].Lineage;
                output += "\t";
                for (j=0;j< NodesList[i].Total_Nuc_mutation.Count;j++)
                output += NodesList[i].Total_Nuc_mutation[j]+",";
                write.WriteLine(output);
            }
            write.Close();
            return;
        }
        static void CalCHNFreq()//计算每个节点的中国序列子代数
        {
            int i, j, k = 0;
            for(i=NodesList.Count-1;i>=0;i--)
            {
                if(NodesList[i].IsLeaf == true)//该节点是叶子节点
                {
                    NodesList[i].Num_AllDescendants = 1;
                    if (NodesList[i].Location.ToUpper() == "CHINA")
                        NodesList[i].Num_CHNDescendants = 1;
                    if (NodesList[i].NewCHN)
                    {
                        NodesList[i].Num_CHNNewDescendants = 1;
                        k++;
                    }
                }
                //累加上儿子节点的数据
                for(j=0;j<NodesList[i].ChildernIndList.Count;j++)
                {
                    NodesList[i].Num_AllDescendants += NodesList[NodesList[i].ChildernIndList[j]].Num_AllDescendants;
                    NodesList[i].Num_CHNDescendants += NodesList[NodesList[i].ChildernIndList[j]].Num_CHNDescendants;
                    NodesList[i].Num_CHNNewDescendants += NodesList[NodesList[i].ChildernIndList[j]].Num_CHNNewDescendants;
                }
            }
            return;
        }
        static string FindAllDescendantNodes(int nodeIndx)//给定节点号，返回所有子孙节点的号 一个字符串 以，分割
        {
            string inds = "" + Convert.ToString(nodeIndx);
            int i, j, k;
            for(j = 0; j < NodesList[nodeIndx].ChildernIndList.Count; j++)
            {
                inds += "," + FindAllDescendantNodes(NodesList[nodeIndx].ChildernIndList[j]);
            }
            return inds;
        }
        static string FindCandiateCHNSubtree(int nodeIndx)//寻找比例大于80%的候选节点，遍历每个节点，如果遇到符合要求的则放弃遍历其子代
        {
            int i, j, k;
            string tmps = "";
            if ((double)NodesList[nodeIndx].Num_CHNNewDescendants / NodesList[nodeIndx].Num_AllDescendants >= 0.8 || (double)NodesList[nodeIndx].Num_CHNDescendants / NodesList[nodeIndx].Num_AllDescendants >= 0.8)
            {
                if (!NodesList[nodeIndx].IsLeaf)
                {
                    tmps += "," + Convert.ToString(nodeIndx);
                }
            }
            else
            {
                for (j = 0; j < NodesList[nodeIndx].ChildernIndList.Count; j++)
                {
                    tmps += FindCandiateCHNSubtree(NodesList[nodeIndx].ChildernIndList[j]);
                }
            }
            return tmps;
        }
        static string LineageSum(string[] nodes)//统计列表中节点的lineage分布
        {
            List<string> LineList = new List<string>();
            List<int> LineCount = new List<int>();
            int i, j, k;
            for(i=0;i<nodes.Count();i++)
            {
                if(NodesList[Convert.ToInt32(nodes[i])].IsLeaf)
                {
                    if (LineList.Contains(NodesList[Convert.ToInt32(nodes[i])].Lineage))
                        LineCount[LineList.IndexOf(NodesList[Convert.ToInt32(nodes[i])].Lineage)]++;
                    else
                    {
                        LineList.Add(NodesList[Convert.ToInt32(nodes[i])].Lineage);
                        LineCount.Add(1);
                    }
                }
            }
            //sort
            int tmpi;string tmps;
            for(i=0;i<LineList.Count();i++)
                for(j=i+1;j<LineList.Count();j++)
                    if(LineCount[i]<LineCount[j])
                    {
                        tmpi = LineCount[i];LineCount[i] = LineCount[j];LineCount[j] = tmpi;
                        tmps = LineList[i];LineList[i] = LineList[j];LineList[j] = tmps;
                    }
            string output = LineList[0];
            for (i = 1; i < LineList.Count(); i++)
                output += "," + LineList[i];
            output += "\t" + Convert.ToString(LineCount[0]);
            for (i = 1; i < LineList.Count(); i++)
                output += "," + Convert.ToString(LineCount[i]);
            
            return output;
        }
        static string CollectionDateSum(string[] nodes)//统计列表中节点的采样时间分布 只有最早和最晚
        {
            List<string> collectDate = new List<string>();
            int i, j, k;
            for (i = 0; i < nodes.Count(); i++)
            {
                if (NodesList[Convert.ToInt32(nodes[i])].IsLeaf && NodesList[Convert.ToInt32(nodes[i])].CollectionDate!=null)
                {
                    collectDate.Add(NodesList[Convert.ToInt32(nodes[i])].CollectionDate);
                }
            }
            //sort
            collectDate.Sort();
            for (i = 0; i < collectDate.Count; i++)
                if (collectDate[i][0] == '2') break;
            return collectDate[i]+"\t"+collectDate[collectDate.Count() - 1];
        }
        static void CollectionDateSum2(string file, string[] nodes)//统计列表中节点的采样时间分布并输出
        {
            StreamWriter write = new StreamWriter(file);
            List<string> collectDate = new List<string>();
            List<int> count = new List<int>();
            int i, j, k;
            for (i = 0; i < nodes.Count(); i++)
            {
                if (NodesList[Convert.ToInt32(nodes[i])].IsLeaf && NodesList[Convert.ToInt32(nodes[i])].Location == "China")
                {
                    if (collectDate.Contains(NodesList[Convert.ToInt32(nodes[i])].CollectionDate))
                        count[collectDate.IndexOf(NodesList[Convert.ToInt32(nodes[i])].CollectionDate)]++;
                    else
                    {
                        collectDate.Add(NodesList[Convert.ToInt32(nodes[i])].CollectionDate);
                        count.Add(1);
                    }
                }
            }
            //sort

            //write
            write.WriteLine("Date\tCount");
            for (i = 0; i < count.Count(); i++)
                write.WriteLine(collectDate[i] + "\t" + Convert.ToString(count[i]));
            write.Close();
            return;
        }
        static void ITolAnnotation(string fold,int target,string[] nodes)
        {
            StreamWriter writenode = new StreamWriter(fold + "/CHNSubtree/" + NodesList[target].name + ".itol.nodeshape.txt");
            StreamWriter writestrip = new StreamWriter(fold + "/CHNSubtree/" + NodesList[target].name + ".itol.strip.txt");
            int i, j, k;
            string output;

            //外圈
            writestrip.WriteLine("DATASET_COLORSTRIP");
            writestrip.WriteLine("SEPARATOR TAB");
            writestrip.WriteLine("BORDER_WIDTH\t0.5");
            writestrip.WriteLine("COLOR\t#bebada");
            writestrip.WriteLine("DATASET_LABEL\tLineage");
            List<string> lineageList = new List<string>();
            for (i = 0; i < nodes.Count(); i++)
                if (NodesList[Convert.ToInt32(nodes[i])].IsLeaf)
                    if (!lineageList.Contains(NodesList[Convert.ToInt32(nodes[i])].Lineage))
                        lineageList.Add(NodesList[Convert.ToInt32(nodes[i])].Lineage);
            string[] ColorsStrip = ColorPicker(lineageList.Count()).Split(',');
            
            output = "LEGEND_COLORS";
            for (i = 0; i < ColorsStrip.Count(); i++)
                output += "\t" + ColorsStrip[i];
            writestrip.WriteLine(output);
            
            output = "LEGEND_LABELS";
            for (i = 0; i < lineageList.Count(); i++)
                output += "\t" + lineageList[i];
            writestrip.WriteLine(output);
            
            output = "LEGEND_SHAPES";
            for (i = 0; i < lineageList.Count(); i++)
                output += "\t1";
            writestrip.WriteLine(output);
            writestrip.WriteLine("LEGEND_TITLE\tLineage");
            writestrip.WriteLine("MARGIN\t5");
            writestrip.WriteLine("STRIP_WIDTH\t35");
            writestrip.WriteLine("DATA");

            for(i=0;i<nodes.Count();i++)
            {
                string name = NodesList[Convert.ToInt32(nodes[i])].name;
                if (NodesList[Convert.ToInt32(nodes[i])].IsLeaf)
                    writestrip.WriteLine(name + "\t" + ColorsStrip[lineageList.IndexOf(NodesList[Convert.ToInt32(nodes[i])].Lineage)] +"\t" + NodesList[Convert.ToInt32(nodes[i])].Lineage);
            }
            writestrip.Close();

            //节点
            writenode.WriteLine("DATASET_SYMBOL");
            writenode.WriteLine("SEPARATOR\tTAB");
            writenode.WriteLine("COLOR\t#6a3d9a");
            writenode.WriteLine("DATASET_LABEL\tCollection Location");
            List<string> LocationList = new List<string>();
            for (i = 0; i < nodes.Count(); i++)
                if (NodesList[Convert.ToInt32(nodes[i])].IsLeaf)
                    if (!LocationList.Contains(NodesList[Convert.ToInt32(nodes[i])].Location))
                        LocationList.Add(NodesList[Convert.ToInt32(nodes[i])].Location);
            string[] ColorsNodes = ColorPicker(LocationList.Count()).Split(',');
            
            output = "LEGEND_COLORS";
            for (i = 0; i < ColorsNodes.Count(); i++)
                output += "\t" + ColorsNodes[i];
            writenode.WriteLine(output);

            output = "LEGEND_LABELS";
            for (i = 0; i < LocationList.Count(); i++)
                output += "\t" + LocationList[i];
            writenode.WriteLine(output);

            output = "LEGEND_SHAPES";
            for (i = 0; i < LocationList.Count(); i++)
                output += "\t2";
            writenode.WriteLine(output);

            writenode.WriteLine("LEGEND_TITLE\tCollection Location");
            writenode.WriteLine("MAXIMUM_SIZE\t10");
            writenode.WriteLine("DATA");

            for (i = 0; i < nodes.Count(); i++)
            {
                string name = NodesList[Convert.ToInt32(nodes[i])].name;
                if (NodesList[Convert.ToInt32(nodes[i])].IsLeaf)
                    writenode.WriteLine(name + "\t2\t10\t" + ColorsNodes[LocationList.IndexOf(NodesList[Convert.ToInt32(nodes[i])].Location)] + "\t1\t1");
            }
            writenode.Close();
            return;
        }
        static void WriteMutationIncidence(string head,string file, string[] nodes, bool ifchina)//计算突变发生率
        {
            StreamWriter write = new StreamWriter(file);
            StreamWriter write_Date = new StreamWriter(file + ".Date");
            write.WriteLine("Node\tTotalDesNum\tCHNDesNum\tNewCHNDesNum\tMutation\tAAMut\tIncidence\tLongestGeneration\tSequenceCount");
            write_Date.WriteLine("Node\tMutation\tAAMut\tDate");
            int i, j, k;
            List<string> MutationList = new List<string>();
            List<string> AAMutList = new List<string>();
            List<int> MutationIncidence = new List<int>();
            List<int> MutationNodeMutNumber = new List<int>();
            List<int> MutationGenerations = new List<int>();
            List<int> MutationSeqCount = new List<int>();

            //先统计背景突变
            Dictionary<string, int> BackgroundMutation = new Dictionary<string, int>();

            //找到中国发生的突变 统计发生次数
            for(i=0;i<nodes.Length;i++)
            {
                int target = Convert.ToInt32(nodes[i]);

                //先统计背景突变
                for (j = 0; j < NodesList[target].Total_Nuc_mutation.Count; j++)
                    if (BackgroundMutation.ContainsKey(NodesList[target].Total_Nuc_mutation[j]))
                        BackgroundMutation[NodesList[target].Total_Nuc_mutation[j]]++;
                    else
                        BackgroundMutation.Add(NodesList[target].Total_Nuc_mutation[j], 1);

                if(!ifchina || (NodesList[target].Location == "China" && NodesList[NodesList[target].FatherInD].Location == "China"))
                {
                    for(j=0;j<NodesList[target].Nuc_mutation.Count();j++)
                    {
                        if (referenceGenome[Convert.ToInt32(NodesList[target].Nuc_mutation[j].Substring(1, NodesList[target].Nuc_mutation[j].Length - 2))] != NodesList[target].Nuc_mutation[j][NodesList[target].Nuc_mutation[j].Length - 1])//不是回复突变
                        {
                            if (MutationList.Contains(NodesList[target].Nuc_mutation[j]))
                            {
                                MutationIncidence[MutationList.IndexOf(NodesList[target].Nuc_mutation[j])]++;
                                if (MutationNodeMutNumber[MutationList.IndexOf(NodesList[target].Nuc_mutation[j])] > NodesList[target].Total_Nuc_mutation.Count())
                                    MutationNodeMutNumber[MutationList.IndexOf(NodesList[target].Nuc_mutation[j])] = NodesList[target].Total_Nuc_mutation.Count();
                                write_Date.WriteLine(head + "\t" + NodesList[target].Nuc_mutation[j] + "\t" + AAMutList[MutationList.IndexOf(NodesList[target].Nuc_mutation[j])] + "\t" + NodesList[target].CollectionDate);
                            }
                            else
                            {
                                MutationList.Add(NodesList[target].Nuc_mutation[j]);
                                MutationIncidence.Add(1);
                                MutationSeqCount.Add(0);
                                MutationGenerations.Add(0);
                                MutationNodeMutNumber.Add(NodesList[target].Total_Nuc_mutation.Count());

                                AAMutList.Add(MutationCombineAnnotation(NodesList[target].Nuc_mutation[j], NodesList[target].Total_Nuc_mutation));
                                write_Date.WriteLine(head + "\t" + NodesList[target].Nuc_mutation[j] + "\t" + AAMutList[MutationList.IndexOf(NodesList[target].Nuc_mutation[j])] + "\t" + NodesList[target].CollectionDate);
                            }
                        }
                    }
                }
            }

            //统计含有突变的序列数量 和 突变最长传递代数
            for (i=0;i<nodes.Length;i++)
            {
                int target = Convert.ToInt32(nodes[i]);
                if (NodesList[target].IsLeaf && NodesList[target].Location == "China")
                { 
                    for (j = 0; j < NodesList[target].Total_Nuc_mutation.Count(); j++)
                    {
                        if (MutationList.Contains(NodesList[target].Total_Nuc_mutation[j]))
                        {
                            MutationSeqCount[MutationList.IndexOf(NodesList[target].Total_Nuc_mutation[j])]++;
                            if (NodesList[target].Total_Nuc_mutation.Count() - MutationNodeMutNumber[MutationList.IndexOf(NodesList[target].Total_Nuc_mutation[j])] > MutationGenerations[MutationList.IndexOf(NodesList[target].Total_Nuc_mutation[j])])
                                MutationGenerations[MutationList.IndexOf(NodesList[target].Total_Nuc_mutation[j])] = NodesList[target].Total_Nuc_mutation.Count() - MutationNodeMutNumber[MutationList.IndexOf(NodesList[target].Total_Nuc_mutation[j])];
                        }
                    }
                }
            }

            //输出
            for(i=0;i<MutationList.Count();i++)
            {
                //忽略背景突变
                if(BackgroundMutation[MutationList[i]]/nodes.Length < 0.95)
                    write.WriteLine(head + "\t" + MutationList[i] + "\t" + AAMutList[i] + "\t" + Convert.ToString(MutationIncidence[i]) + "\t" + Convert.ToString(MutationGenerations[i]) + "\t" + MutationSeqCount[i]);
            }

            write.Close();
            write_Date.Close();
            return;
        }
        static string MutationCombineAnnotation(string mut,List<string> background)
        {
            int i, j, k;
            char[] targetGenome = referenceGenome.ToCharArray();
            background.Add(mut);
            for(i=0;i<background.Count();i++)
            {
                char mutnuc = background[i].Substring(background[i].Length - 1, 1)[0];
                int mutpos = Convert.ToInt32(background[i].Substring(1, background[i].Length - 2));
                targetGenome[mutpos] = (char)mutnuc;
            }

            int pos = Convert.ToInt32(mut.Substring(1, mut.Length - 2));
            string geneName = "TBD";
            int start = 0;
            int end = 0;
            for(i=0;i<GeneList.Count;i++)
            {
                string[] line1 = GeneList[i].Split('\t');
                if(Convert.ToInt32(line1[1]) <= pos && pos <= Convert.ToInt32(line1[2]))
                {
                    geneName = line1[0];
                    start = Convert.ToInt32(line1[1]);
                    end = Convert.ToInt32(line1[2]);
                    break;
                }
            }

            if (i >= GeneList.Count)
                return "intergenic";

            if((pos - start)%3==0)//密码子第一位
            {
                string refcodon = "" + referenceGenome[pos];
                refcodon += referenceGenome[pos + 1];
                refcodon += referenceGenome[pos + 2];
                string mutcodon = new string(targetGenome, pos, 3);
                return geneName + "_" + mimazi[refcodon] + Convert.ToString((pos - start) / 3 + 1) + mimazi[mutcodon];
            }
            if ((pos - start) % 3 == 1)//密码子第2位
            {
                string refcodon = "" + referenceGenome[pos - 1];
                refcodon += referenceGenome[pos];
                refcodon += referenceGenome[pos + 1];
                string mutcodon = new string(targetGenome, pos - 1, 3);
                return geneName + "_" + mimazi[refcodon] + Convert.ToString((pos - start) / 3 + 1) + mimazi[mutcodon];
            }
            if ((pos - start) % 3 == 2)//密码子第3位
            {
                string refcodon = "" + referenceGenome[pos - 2];
                refcodon += referenceGenome[pos - 1];
                refcodon += referenceGenome[pos];
                string mutcodon = new string(targetGenome, pos - 2, 3);
                return geneName + "_" + mimazi[refcodon] + Convert.ToString((pos - start) / 3 + 1) + mimazi[mutcodon];
            }
            return "ERROR";
        }
        static void WriteSeqMutAndDate(string file, string[] nodes)
        {
            StreamWriter write = new StreamWriter(file);
            write.WriteLine("Seq\tDate\tMut\tMutCount\tSpikeAAMut\tSpikeAAMutCount");
            int i, j, k;
            for (i = 0; i < nodes.Count(); i++)
            {
                //if (NodesList[Convert.ToInt32(nodes[i])].IsLeaf && NodesList[Convert.ToInt32(nodes[i])].Location != "China")
                if (NodesList[Convert.ToInt32(nodes[i])].IsLeaf && NodesList[Convert.ToInt32(nodes[i])].Location == "China" && NodesList[Convert.ToInt32(nodes[i])].NewCHN)
                {
                    int target = Convert.ToInt32(nodes[i]);
                    string output = NodesList[target].name;
                    output += "\t" + NodesList[target].CollectionDate;
                    output += "\t" + string.Join(",", NodesList[target].Total_Nuc_mutation);
                    output += "\t" + Convert.ToString(NodesList[target].Total_Nuc_mutation.Count());

                    /*string RBD = "";
                    int RBDCount = 0;
                    for (k = 0; k < NodesList[target].Total_Nuc_mutation.Count(); k++)
                    {
                        if (NodesList[target].Total_Nuc_mutation.Count > 200)
                            Console.WriteLine("?");
                        string aa = MutationCombineAnnotation(NodesList[target].Total_Nuc_mutation[k], NodesList[target].Total_Nuc_mutation);
                        if (aa.Contains("S_"))
                        {
                            RBD += "," + aa;
                            RBDCount++;
                        }
                    }

                    output += "\t" + RBD;
                    output += "\t" + Convert.ToString(RBDCount);*/
                    write.WriteLine(output);
                }
            }
            write.Close();
            return;
        }
        static void FindCHNSubtree()//找usher树上所有的符合要求的中国子树
        {
            int i, j, k;
            string[] CandidateList = FindCandiateCHNSubtree(0).Split(',');
            StreamWriter write = new StreamWriter(JFile + ".CHNClade");
            string output = "Name\tInd\tFatherNode\tSameSeqNum\tTotalDesNum\tCHNDesNum\tNewCHNDesNum\tNewCHNFrac\tLineage\tLineageDist\tEarliestCollectionDate\tLatestCollectionDate\tMutations\n";
            write.Write(output);
            for (i=0;i<CandidateList.Count();i++)
            {
                if (CandidateList[i] != "")
                {
                    int target = Convert.ToInt32(CandidateList[i]);
                    output = NodesList[target].name + "\t";
                    output += CandidateList[i] + "\t";
                    output += Convert.ToString(NodesList[NodesList[target].FatherInD].name) + "\t";
                    output += Convert.ToString(NodesList[target].SameSequenceNumber) + "\t";
                    output += Convert.ToString(NodesList[target].Num_AllDescendants) + "\t";
                    output += Convert.ToString(NodesList[target].Num_CHNDescendants) + "\t";
                    output += Convert.ToString(NodesList[target].Num_CHNNewDescendants) + "\t";
                    output += Convert.ToString((double)NodesList[target].Num_CHNNewDescendants / NodesList[target].Num_AllDescendants) + "\t";

                    //找子代情况
                    string[] Deslist = FindAllDescendantNodes(target).Split(',');
                    output += LineageSum(Deslist) + "\t";
                    output += CollectionDateSum(Deslist) + "\t";

                    //根节点突变
                    output += string.Join(",", NodesList[target].Total_Nuc_mutation.ToArray());

                    //输出iTol注释文件
                    if(NodesList[target].Num_AllDescendants >= 1000)
                    {
                        CollectionDateSum2("M://China220701_230531/CHNSubtree/" + NodesList[target].name + ".date.txt", Deslist);
                        ITolAnnotation("M://China220701_230531", target, Deslist);
                        //输出序列时间 突变
                        WriteSeqMutAndDate("M://China220701_230531/CHNSubtree/" + NodesList[target].name + ".MutDate.txt", Deslist);
                    }

                    //输出突变发生率文件 相较于根节点的回复突变不予考虑
                    if (NodesList[target].Num_AllDescendants >= 1000)
                    {
                        string head = NodesList[target].name + "\t" + Convert.ToString(NodesList[target].Num_AllDescendants) + "\t" + Convert.ToString(NodesList[target].Num_CHNDescendants) + "\t" + Convert.ToString(NodesList[target].Num_CHNNewDescendants);
                        
                        WriteMutationIncidence(head, "M://China220701_230531/CHNSubtree/MutIncidence/" + NodesList[target].name + ".MutIncidence", Deslist, true);
                    }
                    write.Write(output + "\n");
                }
            }
            write.Close();
        }
        static void FindSpecificSubtree()//输出指定节点下的子树的信息
        {
            int i, j, k;
            List<string> SubtreeNameList = new List<string>();
            //SubtreeNameList.Add("node_1084141");
            //SubtreeNameList.Add("node_1035126");
            //SubtreeNameList.Add("node_1093818");
            //SubtreeNameList.Add("node_1024330");
            
            string tmps = "";
            for(i=0;i<NodesList.Count();i++)
            {
                if (SubtreeNameList.Contains(NodesList[i].name))
                    tmps += Convert.ToString(i) + ",";
            }
            string[] CandidateList = tmps.Split(',');
            for (i = 0; i < CandidateList.Count(); i++)
            {
                if (CandidateList[i] != "")
                {
                    int target = Convert.ToInt32(CandidateList[i]);
                    
                    //找子代情况
                    string[] Deslist = FindAllDescendantNodes(target).Split(',');

                    //输出iTol注释文件
                    if (NodesList[target].Num_AllDescendants >= 100)
                    {
                        CollectionDateSum2("M://China220701_230531/CHNSubtree/" + NodesList[target].name + ".date.txt", Deslist);
                        ITolAnnotation("M://China220701_230531", target, Deslist);
                        //输出序列时间 突变
                        WriteSeqMutAndDate("M://China220701_230531/CHNSubtree/" + NodesList[target].name + ".MutDate.txt", Deslist);
                    }

                    //输出突变发生率文件
                    if (NodesList[target].Num_AllDescendants >= 100)
                    {
                        string head = NodesList[target].name + "\t" + Convert.ToString(NodesList[target].Num_AllDescendants) + "\t" + Convert.ToString(NodesList[target].Num_CHNDescendants) + "\t" + Convert.ToString(NodesList[target].Num_CHNNewDescendants);
                        WriteMutationIncidence(head, "M://China220701_230531/CHNSubtree/" + NodesList[target].name + ".MutIncidence", Deslist, true);
                    }
                }
            }
            return;
        }
        static void FindTMRCAForSpecificSubtree()//输出给定子树们的最近祖先节点
        {
            int i, j, k;
            List<string> SubtreeNameList = new List<string>();
            List<int> SubtreeIdList = new List<int>();
            SubtreeNameList.Add("node_1376295");
            SubtreeNameList.Add("node_1352180");
            SubtreeNameList.Add("node_1344263");
            for (i = 0; i < NodesList.Count; i++)
                if (SubtreeNameList.Contains(NodesList[i].name))
                    SubtreeIdList.Add(i);

            List<int> nodePathway = new List<int>();
            j = NodesList[SubtreeIdList[0]].FatherInD;
            while(j!=-1)
            {
                nodePathway.Add(j);
                j = NodesList[j].FatherInD;
            }
            for(i=1;i<SubtreeIdList.Count;i++)
            {
                j = NodesList[SubtreeIdList[i]].FatherInD;
                while (j != -1)
                {
                    if(nodePathway.Contains(j))
                    {
                        nodePathway.RemoveRange(0, nodePathway.IndexOf(j));
                        break;
                    }
                    j = NodesList[j].FatherInD;
                }
            }
            StreamWriter write = new StreamWriter(Workfold + "/CHNSubtree/TheRootOfSpecificSubtrees.tsv");
            write.Write(Convert.ToString(nodePathway[0]) + "\t" + NodesList[nodePathway[0]].name + "\n");
            write.Close();
        }
        static void WriteMutationIncidence_Lineage()//计算同一lineage下非中国序列的突变发生率
        {
            StreamReader read = new StreamReader("M://China220701_230531/CHNSubtree/MutIncidence/LineageList.tsv");
            string line = read.ReadLine();
            List<string> lineageList = new List<string>();
            while(line!=null)
            {
                lineageList.Add(line);
                line = read.ReadLine();
            }
            read.Close();

            int i, j, k=0;
            for(i=0;i<lineageList.Count;i++)
            {
                if (!File.Exists("M://China220701_230531/CHNSubtree/MutIncidence/" + lineageList[i] + ".MutIncidence"))
                {
                    //找到这个lineage下的节点
                    string nodestring = "";
                    for (j = 0; j < NodesList.Count; j++)
                        if (NodesList[j].Lineage == lineageList[i] && NodesList[j].Location != "China")
                        {
                            nodestring += Convert.ToString(j) + ",";
                            k++;
                        }
                    nodestring = nodestring.Substring(0, nodestring.Length - 1);
                    string[] Deslist = nodestring.Split(',');

                    //输出突变率
                    string head = lineageList[i] + "\t" + Convert.ToString(k) + "\tNA" + "\tNA";
                    WriteMutationIncidence(head, "M://China220701_230531/CHNSubtree/MutIncidence/" + lineageList[i] + ".MutIncidence", Deslist, false);
                }
            }
            return;
        }
        static void Main(string[] args)
        {
            JFile = "M://China220701_230531/Data/global_assignments.json";
            Workfold = "M://China220701_230531";
            ReadMetadata();
            ReadinAnnoFile(Workfold);
            MutationEventCal(JFile);
            WriteNodeInfo(JFile);

            CalCHNFreq();//Calculate the Chinese sequence subalgebra for each node
            FindCHNSubtree();//Find all Chinese subtrees on the Usher tree that meet the requirements
            WriteMutationIncidence_Lineage();//Calculate the mutation incidence rate of non Chinese sequences under the same lineage
            FindSpecificSubtree();//Output information of subtrees under specified nodes
            FindTMRCAForSpecificSubtree();//Output to the nearest ancestor nodes of the subtrees
            return;
        }
    }
}

