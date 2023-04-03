import pandas as pd


class ExperimentalData:

    def __init__(self):
        self.data = pd.DataFrame()
        self.time = []
        self.hb = []
        self.hm = []
        self.hz = []

    def _convert_to_dataframe(self):
        # Add Haemoglobin to dataframe
        df_temp = pd.DataFrame(self.hb, columns=['Hb', 'Hb:SEM'], index=self.time)
        self.data = pd.concat([self.data, df_temp])

        # Add free haem to dataframe
        df_temp = pd.DataFrame(self.hm, columns=['Hm', 'Hm:SEM'], index=self.time)
        self.data = self.data.merge(df_temp, left_index=True, right_index=True)

        # Add haemozoin to dataframe
        df_temp = pd.DataFrame(self.hz, columns=['Hz', 'Hz:SEM'], index=self.time)
        self.data = self.data.merge(df_temp, left_index=True, right_index=True)

    def no_drug_nf54(self):
        self.time = [21, 24, 27, 30, 33, 36, 39, 41, 44, 47]  # In hours

        # Haemoglobin: fg/cell, SEM
        self.hb = [[26.5546875, 3.515625], [21.609375, 1.5], [11.6015625, 1.3828125], [4.1015625, 0], [1.4765625, 0],
                   [3.4453125, 0.9140625], [2.296875, 1.3828125], [2.296875, 1.8046875], [1.40625, 0.5390625],
                   [1.6875,	0.9375]]

        # Free haem: fg/cell, SEM
        self.hm = [[1.430397727, 0.183238636], [1.630681818, 0.213068182], [2.393465909, 0.323863636],
                   [2.457386364, 0.136363636], [2.921875, 0.191761364], [3.936079545, 0.289772727],
                   [5.142045455, 0.098011364], [6.066761364, 0.519886364], [5.704545455, 0.234375],
                   [5.815340909, 0.524147727]]

        # Haemozoin: # fg/cell, SEM
        self.hz = [[6.089965398, 2.76816609], [10.51903114,	0], [23.99077278, 3.875432526], [30.08073818, 3.783160323],
                   [35.80161476, 4.290657439], [40.23068051, 4.705882353], [51.30334487, 0], [56.83967705, 6.551326413],
                   [65.88235294, 2.952710496], [69.38869666, 8.581314879]]

        # Convert to dataframe
        self._convert_to_dataframe()

    def no_drug_dd2(self):
        self.time = [20, 23, 26, 29, 32, 35, 38, 41, 44]  # In hours

        # Haemoglobin: fg/cell, SEM
        self.hb = [[1.217, 0.2035], [1.441, 0.3788], [2.107, 0.3346], [1.649, 0.4148], [1.678, 0.141], [2.16, 0.3945],
                   [2.174, 0.7867], [2.459, 0.3721], [1.976, 0.8357]]

        # Free haem: fg/cell, SEM
        self.hm = [[1.48, 0.033], [3.115, 0.9283], [2.082, 0.2012], [2.272, 0.3602], [2.492, 0.1986], [3.878, 0.2985],
                   [3.95, 0.8316], [5.171, 0.8572], [5.83,	1.699]]

        # Haemozoin: # fg/cell, SEM
        self.hz = [[23.05, 8.49], [28.22, 8.72], [31.29, 7.3], [33.76, 4.125], [43.75, 5.025], [54.6, 7.592],
                   [66.76, 12.63], [76.97, 15.91], [97.65, 20.75]]

        # Convert to dataframe
        self._convert_to_dataframe()
