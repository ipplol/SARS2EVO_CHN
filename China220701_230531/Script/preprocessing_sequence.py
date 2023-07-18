from utils import *

def matching_path_pattern(pattern):
    if type(pattern)==type(Path()):
        return glob.glob(str(pattern))
    elif type(pattern)== list:
        filenames = []
        for i in pattern:
            filenames += glob.glob(str(i))
        return filenames

def merge_fasta(path_in, path_out):
    fasta_list = []
    for file in matching_path_pattern(path_in):
        print(file)
        with open(file, 'r') as f:
            fasta = f.read()
        fasta_list.append(fasta)

    fasta_cat = '\n'.join(fasta_list)
    with open(path_out, 'w') as f:
        f.write(fasta_cat)

def merge_metatable(path_in, path_out):
    meta_list = []
    for file in matching_path_pattern(path_in):
        print(file)
        meta = pd.read_csv(file, sep='\t')
        meta_list.append(meta)

    meta_cat = pd.concat(meta_list, axis=0, ignore_index=True)
    meta_cat.to_csv(path_out, sep='\t', index=False)
    return meta_cat

def merge_metatable_selected(path_in, path_out, col):
    meta_list = []
    for file in matching_path_pattern(path_in):
        print(file)
        meta = pd.read_csv(file, sep='\t')
        meta = meta[col]
        meta_list.append(meta)

    meta_cat = pd.concat(meta_list, axis=0, ignore_index=True)
    meta_cat.to_csv(path_out, sep='\t', index=False)

def filter_columns(path_in, path_out):
    df_list = matching_path_pattern(path_in)
    name_list = [i.split('/')[-1] for i in df_list]
    df_list = [pd.read_csv(i, sep='\t',
                           usecols=['strain', 'genbank_accession', 'date',
                                    'country', 'pangolin_lineage'])
               for i in df_list]
    for name, df in zip(name_list, df_list):
        df.to_csv(path_out / f'filter_{name}', sep='\t', index=False)

def rename_fasta(path_in, path_out, name_dic):
    records = SeqIO.parse(path_in, 'fasta')
    records_renamed = []
    for i, record in enumerate(records):
        id = record.id
        description = record.description
        seq = record.seq
        if id in name_dic.keys():
            records_renamed.append(SeqRecord(seq, id=name_dic[id], description=''))
        else:
            # print(id)
            pass

    SeqIO.write(records_renamed, path_out, "fasta")

    return records_renamed

def filter_fasta(path_in, path_out, filter):
    records = SeqIO.parse(path_in, 'fasta')
    records_filtered = []
    id_list = []
    for i, record in enumerate(records):
        id = record.id
        description = record.description
        seq = record.seq
        if id in filter and id not in id_list:
            records_filtered.append(SeqRecord(seq, id=id, description=''))
            id_list.append(id)
        else:
            pass
            #print(id)

    SeqIO.write(records_filtered, path_out, "fasta")

    return records_filtered


def process_ngdc():

    path = Path('/media/data/fuhaoyi/sars2genome/data/china_20230531/')
    path_gisaid_comb = Path(path/'gisaid/combine')
    path_ngdc_raw = Path(path/'ngdc/rawdata')
    path_ngdc_comb = Path(path/'ngdc/combine')
    path_gisaid_ngdc = Path(path/'combine_gisaid_ngdc')

    gisaid_meta = pd.read_csv(path_gisaid_comb/'merge_meta.tsv', sep='\t')
    gisaid_seq = SeqIO.parse(path_gisaid_comb/'merge_seq.fasta', "fasta")
    gisaid_seq = {i.id: i.seq for i in gisaid_seq}


    ngdc_meta = pd.read_csv(path_ngdc_raw / 'metadata_20220701_20230531.tsv', sep='\t')
    # rename meta columns
    # data redundancy elimination
    name = {
        'Virus Strain Name': 'strain',
        'Lineage': 'pangolin_lineage',
        'Sample Collection Date': 'date',
        'Location': 'division',
        'Submitting Lab': 'submitting_lab',
    }
    ngdc_meta = ngdc_meta.rename(columns=name)
    print(ngdc_meta.shape)

    merge_fasta(path_ngdc_raw/'*'/'*.fasta', path_ngdc_comb/'merge_seq.fasta')
    id2name = ngdc_meta.set_index('Accession ID')['strain'].to_dict()
    id2name = {k.split('.')[0]: v for k, v in id2name.items()}
    records = rename_fasta(path_ngdc_comb / 'merge_seq.fasta', path_ngdc_comb / 'merge_seq.fasta', id2name)

    ngdc_seq = SeqIO.parse(path_ngdc_comb/'merge_seq.fasta', "fasta")
    ngdc_seq = {i.id: i.seq for i in ngdc_seq}
    print(len(ngdc_seq))

    # deduplicate
    ngdc_meta = ngdc_meta.drop_duplicates('strain')
    ngdc_meta = ngdc_meta[ngdc_meta['Related ID'].isna()]
    flag = ngdc_meta['date'].apply(lambda x: False if len(x)<8 else True)
    ngdc_meta = ngdc_meta[flag]
    ngdc_meta = ngdc_meta.reset_index(drop=True)

    print(ngdc_meta.shape)
    filter_bool = [True]*ngdc_meta.shape[0]
    gisaid_strain = gisaid_meta['strain'].to_list()
    gisaid_date = gisaid_meta['date'].to_list()
    gisaid_lab = gisaid_meta['submitting_lab'].to_list()
    gisaid_lab = [' '.join(i.split()[0:5]) for i in gisaid_lab]
    identifier_list = zip(gisaid_date, gisaid_lab, gisaid_strain)
    identifier_list = [f'{i[0]}_{i[1]}_{gisaid_seq[i[2]]}' for i in identifier_list]

    for index, row in ngdc_meta.iterrows():
        strain = row['strain']
        date = row['date']
        lab = row['submitting_lab']
        lab = ' '.join(lab.split()[0:5])
        if strain in ngdc_seq.keys():
            identifier = f'{date}_{lab}_{ngdc_seq[strain]}'
            if strain in gisaid_strain:
                filter_bool[index] = False
            if identifier in identifier_list:
                filter_bool[index] = False
        else:
            filter_bool[index] = False

    ngdc_meta = ngdc_meta[filter_bool]
    ngdc_meta = ngdc_meta.reset_index(drop=True)
    print(ngdc_meta.shape)
    ngdc_meta.to_csv(path_ngdc_comb / 'metadata_20220701_20230531_rename_deduplication.tsv', sep='\t', index=False)

    records = filter_fasta(path_ngdc_comb/'merge_seq.fasta', path_ngdc_comb/'merge_seq.fasta', ngdc_meta['strain'].to_list())
    print(len(list(records)))


    # combine gisaid with ngdc
    merge_fasta([path_gisaid_comb/'merge_seq.fasta', path_ngdc_comb/'merge_seq.fasta'], path_gisaid_ngdc/'merge_seq.fasta')
    gisaid_ngdc_meta = merge_metatable([path_gisaid_comb/'merge_meta.tsv', path_ngdc_comb/'metadata_20220701_20230531_rename_deduplication.tsv'],
                    path_gisaid_ngdc/'metadata_20220701_20230531_rename_deduplication.tsv')


    #-------labs---------
    labs = ngdc_meta['submitting_lab'].to_list()
    labs = Counter(labs).items()
    labs = sorted(labs, key=lambda x: x[1], reverse=True)
    labs = pd.DataFrame(labs, columns=['lab', 'count'])
    labs.to_csv(path_ngdc_comb / 'labs.csv', index=False)
    
    labs = gisaid_ngdc_meta['submitting_lab'].to_list()
    labs = Counter(labs).items()
    labs = sorted(labs, key=lambda x: x[1], reverse=True)
    labs = pd.DataFrame(labs, columns=['lab', 'count'])
    labs.to_csv(path_gisaid_ngdc/'labs.csv', index=False)



if __name__ == '__main__':

    # for (root, dirs, files) in os.walk(path_in):
    #     print(root, dirs, files)

    path_in = Path('/media/data/fuhaoyi/sars2genome/data/china_20230531/gisaid/rawdata')
    path_out = Path('/media/data/fuhaoyi/sars2genome/data/china_20230531/gisaid/combine')
    # merge_fasta(path_in/'*'/'*.fasta', path_out/'merge_seq.fasta')
    # merge_metatable(path_in/'*'/'*.tsv', path_out/'merge_meta.tsv')

    process_ngdc()












