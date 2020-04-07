
import statistics
import os
import re
from xml.etree import ElementTree
import xlsxwriter


class ParserPep:

    def __init__(self, input, mode):
        self.file_name = input
        #self.output_file_name = output
        self.dict_peptides = {}
        self.namespace = ""
        self.root = ""
        self._init_root_and_namespace()
        self.c = 0
        self.running_mode = mode

    # there is a namespace before root, it is pain in the ass
    # so this function cut it so other functions in ElementTree would work properly
    def _namespace(self, element):
        m = re.match(r'\{.*\}', element.tag)
        return m.group(0) if m else ''

    # def k_count(self, seq: str):
    #     """
    #     how many lysines(unmarked+heavy) in peptide
    #     """
    #     return seq.count("K")
    # def k_heavy_count(self, seq: str, counter_all):
    #     """
    #     how many heavy lysines (marked K[162]) in peptide
    #     param "counter_all" is the return value of k_count - just for assertion
    #     """
    #     temp = seq.count("K[162]")
    #     assert (temp >= counter_all)
    #     return temp

    def little_test(self, mod: str):
        assert (mod.count("[43]") == 1)

    def _k_peptide_type(self, seq):
        """
        this function is relevant only for "lysine" mode (so far)
        if the peptide have no lysines(K) in it - return "no k"
        if the peptide is heavy == all lysines are [162] - return "heavy"
        id the peptide is light == all lysines are unmarked - return "light"
        if the peptide is mixed (heavy&light) - return "bad kitty"
        """
        k_all = seq.count("K")
        k_heavy = seq.count("K[162]")
        k_light = k_all - k_heavy
        if k_all == 0:
            # no lysines in pep
            return "no k"
        elif k_heavy == k_all:
            assert (k_all > 0)
            assert (k_light == 0)
            return "heavy"
        elif k_light == k_all:
            assert (k_all > 0)
            assert (k_heavy == 0)
            return "light"
        else:
            assert (k_all > 0)
            return "bad kitty"


    # check if there is modification in the n term
    # start with "n[29 |15|18 | 35] then the seq of peptide
    def _modification_n(self, seq, mod):
        #self.c += 1
        if re.match('n\[29\].', mod) != None or re.match('n\[15\].', mod)!= None:
            self.dict_peptides[seq].add_n_light()
        elif re.match('n\[35\].', mod) != None or re.match('n\[18\].', mod)!= None:
            self.dict_peptides[seq].add_n_heavy()


    # error rate is the minimum probability we want
    # if peptide has probability < error_rate we pass it
    def _calc_error_rate(self, err):
        error = self.root.findall('.//' + self.namespace + 'roc_error_data')
        # if we  do not find the filed charge =all we return -1
        error_rate = -1
        for x in error:
            if x.attrib['charge'] == "all":
                error_point = x.findall('.//' + self.namespace + 'error_point')
                for y in error_point:

                    if y.attrib['error'] == str(err):
                        error_rate = y.attrib['min_prob']
                        error_rate = float(error_rate)
                        break
                    else:
                        continue
                else:
                    continue
        return error_rate

    def _init_root_and_namespace(self):
        full_file = os.path.abspath(self.file_name)
        dom = ElementTree.parse(full_file)
        self.root = dom.getroot()
        self.namespace = self._namespace(self.root)



    def parse_dict(self, error_rate):
        error_rate = self._calc_error_rate(error_rate)
        # print (name_space)
        hits = self.root.findall('.//' + self.namespace + 'search_hit')
        #dict_peptides = {}
        counter = 0
        for x in hits:
            probability = x.find('.//' + self.namespace + 'peptideprophet_result')
            if probability != None:
                probability = probability.attrib['probability']
            else:
                continue
            # print(probability)
            probability = float(probability)
            # print(probability)
            if probability < error_rate:
                continue
            if self.running_mode == "default":
                heavy = x.find('.//' + self.namespace + 'xpressratio_result')
                if heavy != None:
                    heavy = heavy.attrib['heavy_area']
                    heavy = float(heavy)
                else:
                    continue
                light = x.find('.//' + self.namespace + 'xpressratio_result')
                if light != None:
                    light = light.attrib['light_area']
                    light = float(light)
                else:
                    continue

                mod = x.find('.//' + self.namespace + 'modification_info')
            elif self.running_mode == "label":
                mod = None
                matched_ions = x.attrib['num_matched_ions']
                label_free_res = x.find('.//' + self.namespace + 'xpresslabelfree_result')
                if label_free_res != None:
                    peak_area = label_free_res.attrib['peak_area']
                    peak_area = float(peak_area)
                    rt_seconds = label_free_res.attrib['peak_intensity_RT_seconds']
                    rt_seconds = float(rt_seconds)
                    peak_intensity = label_free_res.attrib['peak_intensity']
                    peak_intensity = float(peak_intensity)
                else:
                    continue

            # else:
            #     print("shouldnt get here!!")
            #     exit(1)
            if self.running_mode == "lysine":
                mod = x.find('.//' + self.namespace + 'modification_info')
                #print("bloop")
                #print(mod)
                if mod != None:
                    mod = mod.attrib['modified_peptide']
                    #print(mod)
                    k_type = self._k_peptide_type(mod)
                    if k_type == "bad kitty":
                        continue
                    else:
                        heavy = x.find('.//' + self.namespace + 'xpressratio_result')
                        if heavy != None:
                            heavy = heavy.attrib['heavy_area']
                            heavy = float(heavy)
                        elif k_type == "no k":
                            heavy = 0
                        else:
                            continue
                        light = x.find('.//' + self.namespace + 'xpressratio_result')
                        if light != None:
                            light = light.attrib['light_area']
                            light = float(light)
                        elif k_type == "no k":
                            light = 0
                        else:
                            continue
                else:
                    continue
            seq = x.attrib['peptide']
            if seq in self.dict_peptides:
                self.dict_peptides[seq].add_counter()   # in all modes
                if self.running_mode == "lysine":
                    #print("here")
                    self.little_test(str(mod))
                    k_type = self._k_peptide_type(mod)
                    assert (k_type != "bad kitty")
                    self.dict_peptides[seq].add_heavy(heavy)
                    self.dict_peptides[seq].add_light(light)
                    if k_type == "heavy":
                        self.dict_peptides[seq].add_n_heavy()
                    elif k_type == "light":
                        self.dict_peptides[seq].add_n_light()
                    else:
                        assert (k_type == "no k")
                        self.dict_peptides[seq].add_no_k()

                    continue
                if self.running_mode == "default":
                    self.dict_peptides[seq].add_heavy(heavy)
                    self.dict_peptides[seq].add_light(light)
                    if mod != None:
                        mod = mod.attrib['modified_peptide']
                        self._modification_n(seq, mod)
                        #self.dict_peptides[seq].print_peptide()
                        continue
                elif self.running_mode == "label":
                    self.dict_peptides[seq].add_peak_area(peak_area)
                    self.dict_peptides[seq].add_avg_rt_seconds(rt_seconds)
                    self.dict_peptides[seq].add_peak_intensity(peak_intensity)
                    continue
                else:
                    print("shouldnt get here1!!")
                    exit(1)

            pep_type = x.attrib['num_tol_term']
            prot = x.attrib['protein']
            prot_alternative = x.findall('.//' + self.namespace + 'alternative_protein')
            list_alternative = [prot]
            if prot_alternative:
                for y in prot_alternative:
                    curr = y.attrib['protein']
                    if curr not in list_alternative:
                        list_alternative.append(curr)

                # print(list_alternative)
                # print(prot)
            # else:
            #     print("sould be NOne*************************************")
            #     print(prot_alternative)
            #     exit(1)

            start = x.attrib['peptide_prev_aa']
            end = x.attrib['peptide_next_aa']
            if self.running_mode == "default" or self.running_mode == "lysine":
                # -1 for peak_area, peak intensity, rt_seconds, ions
                self.dict_peptides[seq] = Peptide(seq, pep_type, prot, list_alternative, probability, start, end, heavy,
                                                  light, -1, -1, -1, -1, self.running_mode)
                if self.running_mode == "lysine":
                    if(k_type == "no k"):
                        #print(k_type)
                        self.dict_peptides[seq].add_no_k()
            elif self.running_mode == "label":
                # -1 for heavy and light
                self.dict_peptides[seq] = Peptide(seq, pep_type, prot, list_alternative, probability, start, end, -1,
                                                  -1, peak_area,peak_intensity, rt_seconds, matched_ions, self.running_mode)
            else:
                print("shouldnt get here2!!")
                exit(1)
            if self.running_mode == "default":
                if mod != None:
                    #assert (self.running_mode == "default")
                    mod = mod.attrib['modified_peptide']
                    self._modification_n(seq, mod)
                #self.dict_peptides[seq].print_peptide()
            elif self.running_mode == "lysine":
                k_type = self._k_peptide_type(mod)
                assert (k_type != "bad kitty")
                if k_type == "heavy":
                    self.dict_peptides[seq].add_n_heavy()
                elif k_type == "light":
                    self.dict_peptides[seq].add_n_light()


class Peptide:
    """
    parsing pep.xml into dictionary
    the dict has key=seq and value= Peptide class
    """

    def __init__(self, seq, pep_type, prot, alternative, prob, start, end, heavy , light, peak_area, peak_intensity,
                 rt_seconds, ions, mode):
        self.seq: str = seq                   # attrib "peptide" of "search_hit" in xml
        self.pep_type: int = pep_type         # attrib "num_tol_term" of "search_hit in xml
        self.prot: str = prot                 # attrib "protein" of "search_hit" in xml
        self.alternative: list = alternative  # all allternative_protein.attrib(protein)
        self.prob: float = prob               # attrib  "probabilty" of "search_hit\ peptideprophet_result" in xml
        self.start: str = start               # attrib "peptide_prev_aa"  of "search_hit" in xml
        self.end: str = end                   # attrib "pepide_next_aa"  of "search_hit" in xml
        self.counter: int = 1                 # number of occurrences of the same seq in file
        self.n_heavy: int = 0                 # number of peptides with "n[35]" modifications
        self. n_light: int = 0                # number of peptides with "n[29]" modificationa
        self.heavy: float = heavy             # sum of all heavy_area - attrib of search_hit\ xpressratio_result
        self.light: float = light             # sum of all light_area - attrib of search_hit\ xpressratio_result
        self.peak_area: float = peak_area     # sum of all peak area for label-free mode
        self.peak_intensity: float = peak_intensity  # sum of all peak intensity for label-free mode
        self.rt_seconds: float = rt_seconds   # avg og all rt_seconds intensity for label-free mode
        self.ions: int = ions                 # number of matched ions - for label-free mode
        self.ratio: str = "0"                 # light / heavy. if heavy=0 ratio=-1
        self.mode: int = mode                 # is it free- label (==1) or not (==0)
        self.no_k: int = 0                    # used only for uniform n k running mode to sum up the ynmarked peps
        if self.mode == "default" or self.mode == "lysine":
            self._calc_ratio()

    def add_n_heavy(self):
        self.n_heavy += 1

    def add_n_light(self):
        self.n_light += 1

    def add_counter(self):
        self.counter += 1

    def add_heavy(self, new_heavy):
        self.heavy += new_heavy
        # recalcing ratio
        self._calc_ratio()

    def add_light(self, new_light):
        self.light += new_light
        # recalcing ratio
        self._calc_ratio()

    def add_no_k(self):
        self.no_k += 1

    def add_peak_area(self, new_peak):
        self.peak_area += new_peak

    def add_peak_intensity(self, new_intens):
        self.peak_intensity += new_intens

    def add_avg_rt_seconds(self, new_rt):
        sum_without_new = (self.counter - 1) * self.rt_seconds
        self.rt_seconds = (sum_without_new + new_rt) / self.counter

    def _calc_ratio(self):
        if self.heavy == 0 and self.light == 0:  # preventing dev in 0
            self.ratio = "0 dev 0"
        elif self.heavy == 0:
            assert self.light
            self.ratio = "num dev 0"
        else:
            assert self.heavy
            ratio: float = self.light/self.heavy
            str(ratio)
            self.ratio = ratio

    # just for testing
    def print_peptide(self):
        print("seq is " + self.seq)
        #print("type is " + self.pep_type)
        #print("prot is " + self.prot)
        # cannot print srr with none-str type... so i used join to convert list to str
        #print("alternative is " + "".join(self.alternative))
        #print("prob is " + str(self.prob))
        #print("start is " + str(self.start))
        #print("end is " + str(self.end))
        print("counter is " + str(self.counter))
        print("heavy is " + str(self.heavy))
        print("light is " + str(self.light))
        print("ratio is " + str(self.ratio))
        print ("n heavy is " + str(self.n_heavy))
        print("n light is " + str(self.n_light))
        print("**********************")

    def class_to_list(self, mode):
        pep = []
        pep.append(self.seq)
        pep.append(self.pep_type)
        pep.append(self.prot)
        pep.append("".join(self.alternative))
        pep.append(self.prob)
        pep.append(self.start)
        pep.append(self.end)
        pep.append(self.counter)
        if mode == "default" or mode == "lysine":
            pep.append(self.n_heavy)
            pep.append(self.n_light)
            if mode == "lysine":
                pep.append(self.no_k)
            pep.append(self.heavy)
            pep.append(self.light)
            pep.append(self.ratio)
        elif mode == "label":
            pep.append(self.peak_area)
            pep.append(self.peak_intensity)
            pep.append(self.rt_seconds)
            pep.append(self.ions)
        else:
            print("shouldnt get here3!!")
            exit(1)
        return pep


class United:
    """
    could not think about somethimg more convinient or smart
    so i created this class for the merger file
    a new data type to unite all seq from all dict_peptides
    count all appearances in all files
    ratios separately
    """
    def __init__(self, seq, prot, count, num_of_files, mode, start, end):
        self.seq: str = seq
        self.protein: str = prot
        self. sum_count: int = count
        self.start = start
        self.end = end
        self.ratio_dict = {}
        if mode == "default" or mode == "lysine":
            for i in range(num_of_files):
                self.ratio_dict[i] = ""
        elif mode == "label":
            self.count_dict = {}
            for i in range(num_of_files):
                self.count_dict[i] = ""
            self.peak_area_dict = {}
            for i in range(num_of_files):
                self.peak_area_dict[i] = ""
            self.peak_intensity_dict = {}
            for i in range(num_of_files):
                self.peak_intensity_dict[i] = ""
            self.rt_seconds_dict = {}
            for i in range(num_of_files):
                self.rt_seconds_dict[i] = ""
            self.ion_dict = {}
            for i in range(num_of_files):
                self.ion_dict[i] = ""

        self.avg = -1
        self.median = -1
        self.st_deviation = -1
        self.mode = mode

    def add_sum_count(self, count):
        self.sum_count += count

    def add_ratio(self, i, ratio):
        self.ratio_dict[i] = ratio
        self._calc_all()

    def add_all_label_mode(self, i, counter, peak_area, peak_intensity, rt_seconds, ions):
        self.count_dict[i] = counter
        self.peak_area_dict[i] = peak_area
        self.peak_intensity_dict[i] = peak_intensity
        self.rt_seconds_dict[i] = rt_seconds
        self.ion_dict[i] = ions

    def print_united(self):
        print("seq is " +self.seq)
        print("prot is " +self.protein)
        print( "count is " + str(self.sum_count))
        print ("ratio dict is " + str(self.ratio_dict))

    def class_to_list(self, mode):
        pep = []
        pep.append(self.seq)
        pep.append(self.protein)
        pep.append(self.sum_count)
        pep.append(self.start)
        pep.append(self.end)
        if mode == "default" or mode == "lysine":
            for i in range(len(self.ratio_dict)):
                pep.append(self.ratio_dict[i])
            pep.append(self.avg)
            pep.append(self.median)
            pep.append(self.st_deviation)

        elif mode == "label":
            for i in range(len(self.count_dict)):
                pep.append(self.count_dict[i])
            for i in range(len(self.peak_area_dict)):
                pep.append(self.peak_area_dict[i])
            for i in range(len(self.peak_intensity_dict)):
                pep.append(self.peak_intensity_dict[i])
            for i in range(len(self.rt_seconds_dict)):
                pep.append(self.rt_seconds_dict[i])
            for i in range(len(self.ion_dict)):
                pep.append(self.ion_dict[i])
        return pep

    def _dict_to_list_int(self):
        l = []
        for i in range(len(self.ratio_dict)):
            curr = self.ratio_dict[i]
            if curr == "0 dev 0" or curr == "num dev 0" or curr =="":
                continue
            else:
                l.append(self.ratio_dict[i])
        return l

    def _calc_all(self):
        ratio_list = self._dict_to_list_int()
        if len(ratio_list) > 0:
            self.avg = statistics.mean(ratio_list)
            self.median = statistics.median(ratio_list)
            if len(ratio_list) > 1:
                self.st_deviation = statistics.stdev(ratio_list)


class Model:

    def __init__(self):
        self.reader = None

    def validation(self, files, n):
        # n= number of files
        if n == 0:
            # user did not choose files
            return "Exit"
        for i in range(n):
            if not re.search(r"\.xml$", files[i]):
                # if not all files are XML files
                return "Invalid"
        else:
            return "Valid"

    def xlsx_create(self, output_file_name, dict_peptides, header, mode):
        """"
        export the data to excel sheet
        """
        workbook = xlsxwriter.Workbook(output_file_name + '.xlsx')
        sheet = workbook.add_worksheet()
        headers = header
        headers_len = len(headers)
        for x in range(headers_len):
            sheet.write(0, x, headers[x])
        # col and raw are just for the for iterations below
        row = 1
        for pep in dict_peptides:
            ll = dict_peptides[pep].class_to_list(mode)
            if headers_len == 17:
                print("headers len is " + str(headers_len))
                print("ll len is " + str(len(ll)))
                print(headers)

            assert (headers_len == len(ll))
            for x in range(headers_len):

                sheet.write(row, x, ll[x])
            row += 1
        workbook.close()  # creating the output.xlsx

    def header_unite_create(self, num_of_files, mode):
        headers = ["seq", "protein", "counter all", "start", "end"]
        if mode == "default" or mode == "lysine":
            #assert (mode == 0)
            for i in range(num_of_files):
                headers.append("ratio" + str(i + 1))
            headers.append("mean")
            headers.append("median")
            headers.append("st dev")
        elif mode == "label":
            #assert (mode == 1)
            for i in range(num_of_files):
                headers.append("count in file " + str(i + 1))
            for i in range(num_of_files):
                headers.append("peak area" + str(i + 1))
            for i in range(num_of_files):
                headers.append("peak intensity " + str(i + 1))
            for i in range(num_of_files):
                headers.append("rt seconds " + str(i + 1))
            for i in range(num_of_files):
                headers.append("matched ions " + str(i + 1))

        return headers

    def file_parse(self, file_name, output_name, error_rate, mode):
        f = ParserPep(file_name, mode)
        f.parse_dict(error_rate)
        dict_f = f.dict_peptides
        if mode == "lysine":
            headers = ["seq", "type(semi/full)", "protein", "alternative protein", "probability", "start", "end", "counter",
                       "k mod heavy", "k mod light", "no k", "sum heavy", "sum light", "ratio"]
            self.xlsx_create(output_name, dict_f, headers, mode)
        if mode == "default":
            headers = ["seq", "type(semi/full)", "protein", "alternative protein", "probability", "start", "end", "counter",
                       "n mod heavy", "n mod light", "sum heavy", "sum light", "ratio"]
            self.xlsx_create(output_name, dict_f, headers, mode)
        elif mode == "label":
            headers = ["seq", "type(semi/full)", "protein", "alternative protein", "probability", "start", "end",
                       "counter", "sum peak area", "sum peak intensity", "avg rt seconds", "matched ions"]
            self.xlsx_create(output_name, dict_f, headers, mode)

        return dict_f



