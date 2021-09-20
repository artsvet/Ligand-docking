from proteins import Multiprocess, Protein
import ligands
from ligands import Mol, Lipid, Lea, Pl
from vinaDock import Docker


def KD(dG, T = 310):
    return 1 / np.e ** (-dG / (R * T))


DEFAULT_LEAS = ['LEA8', 'LEA16', 'LEA18Z9', 'LEA20Z581114']


class TargetsParse:
    """
    Parses Series to set up dependencies for ligand
    docking. Contains Protein target attributes,
    specifies ligands that will be docked using
    reference lipid and pl tables. Will serve as
    primary interface for running dock and
    collecting outputs.

    todo: abstract dock into class that takes targetsParse
    """

    def __init__(self, dock_series, lipid_table):
        self.lipid_table = lipid_table
        self.id = dock_series['id']
        self.pdb = Protein(dock_series['pdb'])
        self.name = dock_series['name']
        self.type = dock_series['type']
        self.species = dock_series['species']
        self.function = dock_series['function']
        self.lipid_patterns = re.split(', ', dock_series['lipidPatterns'])
        self.lipid_tails = [tails for tails
                            in re.split('; ', dock_series['lipidtails'])
                            ]
        self.selected_names = re.split('  ', dock_series['selectedNames'])
        self.selected_smiles = re.split('  ', dock_series['selectedSmiles'])
        assert len(self.selected_names) == len(self.selected_smiles)
        self.refLigands = [ligand for key in self.lipid_patterns
                           if key != ''
                           for ligand in DEFAULT_LIGANDS[key]
                           ]
        self.box = eval(dock_series['box'])
        self.lea = [LeaFromSeries(series)
                    for index, series in self.lipid_table.iterrows()]
        self.lipids = self.build_lipid_list()
        self.selects = self.build_select_list()

    def build_lipid_list(self):
        lipid_list = []

        if self.lipid_patterns == ['']:
            return []

        for pattern, sidechain in itertools.product(
                self.lipid_patterns, self.lipid_tails
        ):

            if pattern not in PL_PATTERNS:
                continue
            if pattern == 'lipid':
                try:
                    assert ' ' not in sidechain
                except AssertionError:
                    continue
                lipid_list.append(LipidFromString(sidechain))
            else:
                lipid_list.append(PlFromString(pattern, sidechain))

        return lipid_list

    def build_select_list(self):
        select_list = []
        if self.selected_names == ['']:
            return []

        for i, name in enumerate(self.selected_names):
            ligand = Mol(self.selected_smiles[i])
            setattr(ligand, 'name', name)
            select_list.append(ligand)
        return select_list

    def annotate_long(self, longDf):

        speciesMask = 1 if self.species == 'H.sapiens' else 0
        leaMask = 1 if longDf.iloc[0]['Ligand'] in DEFAULT_LEAS else 0
        ligMask = 1 if longDf.iloc[0]['Ligand'][:2] in self.lipid_patterns else 0

        longDf['species'] = np.full(len(longDf.index), speciesMask)
        longDf['defaultLea'] = np.full(len(longDf.index), leaMask)
        longDf['defaultLigand'] = np.full(len(longDf.index), ligMask)

        return longDf

    def dock(self, ligands, run_count=1, cpu=8,
             exhaustiveness=10, write_dir=os.getcwd()):

        long = '{0}{1}_{2}_out.csv'.format(
            write_dir, self.name, self.id)
        trimmed = '{0}{1}_{2}_trimmed.csv'.format(
            write_dir, self.name, self.id)
        short = '{0}{1}_{2}_summary.csv'.format(
            write_dir, self.name, self.id)
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)
        os.chdir(write_dir)

        if os.path.exists(long):
            long_df = pd.read_csv(long)
        else:
            long_df = pd.DataFrame(
                columns=['Date_time', 'Exhaustiveness', 'Run_number',
                         'Receptor', 'Ligand', 'Rank', 'Affinity',
                         'Dist_rmsd_l.b.', 'Dist_rmsd_u.b.', 'species',
                         'defaultLea', 'defaultLigand'
                         ]
            )
            long_df.to_csv(long, index=False, mode='w')

        if os.path.exists(trimmed):
            trim_df = pd.read_csv(trimmed)
        else:
            trim_df = long_df
            trim_df.to_csv(trimmed, index=False, mode='w')

        if os.path.exists(short):
            short_df = pd.read_csv(short, header=None,
                                   names=['', 'dGmean', 'dGsd', 'KDmean', 'KDsd'])
            short_df.to_csv(short, header=False, index=False, mode='w')
        else:
            short_df = pd.DataFrame(
                columns=['', 'dGmean', 'dGsd', 'KDmean', 'KDsd'],
                data=[['', self.id],
                      ['type', self.type],
                      ['species', self.species],
                      ['lipids', ', '.join(self.lipid_patterns)],
                      [],
                      ['', 'dG mean', 'dG sd', 'KD mean', 'KD sd'],
                      ['ref', np.nan, np.nan, np.nan, np.nan]],
            )

        for mol in ligands:

            cpu_lock = time.localtime(time.time()).tm_hour
            if 21 > cpu_lock > 8:
                cpu = 6

            dock_output = self.annotate_long(
                Docker(receptor=self.pdb,
                       ligand=mol,
                       log_path='{0}\\{1}_{2}_log.txt'.format(
                           write_dir, self.id, mol.name),
                       box=self.box, run_count=run_count,
                       exhaustiveness=exhaustiveness,
                       cpu=cpu
                       ).run()
            )

            long_df = long_df.append(dock_output)
            trim = dock_output.loc[dock_output['Rank'] == '1']
            trim_df = trim_df.append(trim)

            dg = np.mean(trim_df.loc[
                             trim_df['Ligand'] == mol.name, 'Affinity'
                         ].astype(float)), \
                 np.std(trim_df.loc[
                            trim_df['Ligand'] == mol.name, 'Affinity'
                        ].astype(float))
            kd = KD(dg[0]), \
                 KD(dg[0]) * (-dg[1] / dg[0])

            short_df.loc[6, 'dGmean'] = (
                np.mean(trim_df.loc[
                            trim_df['defaultLigand'] == 1, 'Affinity'].astype(float))
            )
            short_df.loc[6, 'dGsd'] = (
                np.std(trim_df.loc[
                           trim_df['defaultLigand'] == 1, 'Affinity'].astype(float))
            )
            short_df.loc[6, 'KDmean'] = KD(short_df.loc[6, 'dGmean'])
            short_df.loc[6, 'KDsd'] = short_df.loc[6, 'KDmean'] * \
                                      -(short_df.loc[6, 'dGsd'] / short_df.loc[6, 'dGmean'])
            short_df = short_df.append(
                pd.DataFrame(
                    columns=['', 'dGmean', 'dGsd', 'KDmean', 'KDsd'],
                    data=[[mol.name, dg[0], dg[1], kd[0], kd[1]]]
                ), ignore_index=True
            )

            dock_output.to_csv(long, index=False, header=False, mode='a')
            trim.to_csv(trimmed, index=False, header=False, mode='a')
            short_df.to_csv(short, index=False, header=False, mode='w')

        return long_df, trim_df, short_df