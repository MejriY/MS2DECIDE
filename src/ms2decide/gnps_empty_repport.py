from .ClosestGNPS import _gnps_annotations_download_results
from .ClosestGNPS2 import _gnps2_annotations_download_results
import pandas as pd


def gnps_empty_repport(plateforme, dict_task_id, df,save_path):
    while (plateforme not in ['GNPS' or 'GNPS2']):
        plateforme = input('PLEASE SELECT GNPS or GNPS2')
    if (plateforme == 'GNPS'):
        res = {}
        for peak in dict_task_id:
            s = {i: 0 for i in dict_task_id[peak]}
            for mass in dict_task_id[peak]:
                job = dict_task_id[peak][mass]
                s[mass] = _gnps_annotations_download_results(job)
            res[peak] = s
    else:
        res = {}
        for peak in dict_task_id:
            s = {i: 0 for i in dict_task_id[peak]}
            for mass in dict_task_id[peak]:
                job = dict_task_id[peak][mass]
                s[mass] = _gnps2_annotations_download_results(job)
            res[peak] = s
    x = df[df.inchi_gnps == '*']['ID'].values.tolist()
    url = 'https://gnps.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID={}'
    report = []
    taken = set()
    for ID in list(x):
        # print(ID)
        for peak in [6, 5, 4]:
            for mass in [0.02, 0.1, 10, 25, 50, 100, 250, 500, 1000]:
                t = res[peak][mass]
                l = t[t['#Scan#'] == ID].values.tolist()
                if (len(l) != 0) and (ID not in taken):
                    taken.add(ID)
                    report.append(
                        [ID, peak, mass, l[0][0], url.format(l[0][0])])
    REP = pd.DataFrame(
        report, columns=['ID', 'peak', 'mass', 'CCMSLIB', 'link CCMSLIB'])
    REP.to_csv(save_path, index=False,sep='\t')
