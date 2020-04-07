import tkinter as tk
from tkinter import ttk, DISABLED, NORMAL
import os
from tkinter.filedialog import askopenfilename
from tkinter import messagebox
from PIL import Image, ImageTk
from itertools import count
from xml_parser_model import Model, United
import threading


MAIN_FONT = ("Arial Bold", 10)


class Controller(tk.Tk):
    """
    creates two graphic pages and comunicates with Model- the actual parser
    """
    def __init__(self):
        self.model = Model()            # the actual parser
        self.files: list = list()       # all pep-xml files to parse and unite
        self.swap_dict = dict()         # key: name of file value: 0 if light/heavy and 1 if heavy/light
        self.dict_list = list()         # all parsered files
        self.number_of_files = 0
        tk.Tk.__init__(self)
        #   cat icon misbehaving in linux
        #tk.Tk.iconbitmap(self, default="catIcon.ico")
        self.error_rate_button = tk.StringVar()     # 0/ 0.01 / 0.005 / 0.001
        self.running_options = tk.StringVar()       # "default" "label" or "lysine"
        self.output_entry = ""  # name of united file output
        #self.label_free = tk.IntVar()   # if 1 than consider
        self._initialize_radio_buttons()
        self.title("PEP-XML PARSER")
        self.container = tk.Frame(self)
        self.container.pack(side="top", fill="both", expand=True)
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)
        self.frames: dict = {}
        frame = StartPage(self)
        self.frames[StartPage] = frame
        frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def _print_file_numbers(self):
        for i, f in zip(range(self.number_of_files), self.files):
            print(str(i+1) + "   " + str(f))

    def _initialize_radio_buttons(self):
        #self.label_free.set(0)
        self.running_options.set("-1")
        self.error_rate_button.set("-1")

    def choose_file(self):
        self.files = self.frames[StartPage].browse()
        self.number_of_files = len(self.files)
        validation = self.model.validation(self.files, self.number_of_files)
        if validation == "Exit":
            self.frames[StartPage].change_number_of_files_label(self.number_of_files)
            return
        elif validation == "Invalid":
            self.error_message("Invalid file", "All files must be XMLs")
            return
        else:
            self.frames[StartPage].change_number_of_files_label(self.number_of_files)

    def next(self):
        """
        checking the user input- if everything is ok moving to second page
        """
        if not len(self.files):
            self.error_message("Invalid file", "No file Selected")
            return
        if self.output_entry.get() == "":
            self.error_message("Invalid entry", "Please enter output file name")
            return
        if self.error_rate_button.get() == "-1":
            self.error_message("Invalid choice", "Please choose error rate")
            return
        if self.running_options.get() == "-1":
            self.error_message("Invalid choice", "Please choose running mod")
            return
        self._print_file_numbers()
        # print("output - error- lable-free")
        # print(self.output_entry.get())
        # print(self.error_rate_button.get())
        # print(self.label_free.get())
        self.run()


        # frame = SecondPage(self)
        # self.frames[SecondPage] = frame
        # self.frames[StartPage].grid_remove()    # remove so that tk will initialize the frame geometry
        #                                         # according to the second page
        # frame.grid(row=0, column=0, sticky="nsew")
        # self.show_frame(SecondPage)

    # def back_to_start_page(self):
    #     self.frames[SecondPage].grid_remove()   # remove so that tk will initialize the frame geometry
    #     self.frames[StartPage].grid()
    #     self.show_frame(StartPage)

    def run(self):

        # for d in self.swap_dict:
        #     print(self.swap_dict[d].get())
        for f in self.files:
            #parser files one by one. add to the dict-list after parsing
            #output name is with "_out" ending
            try:
                self.dict_list.append(self.model.file_parse(f, str(f + '_out'), self.error_rate_button.get(),
                                                            self.running_options.get()))
            except:
                # it is kinda early stage to put a try-catch but i could not predict what will go wrong in the parser
                # and i didnt want it to have an interpreter error but a gui one
                # so i catch them all like pokemons
                self.error_message("Invalid file", "Cannot parse\ncould be wrong mode or an open xlsx in use")
                os._exit(1)

        n = len(self.dict_list)
        assert (n == len(self.files))
        print("all ok until united")
        #exit(0)
        merge_dict = {}
        for i, d in zip(range(n), self.dict_list):
            for seq in d:
                if seq not in merge_dict:
                    merge_dict[seq] = United(seq, d[seq].prot, d[seq].counter, n, self.running_options.get(), d[seq].start,
                                             d[seq].end)
                    if self.running_options.get() == "default" or self.running_options.get() == "lysine":
                        merge_dict[seq].add_ratio(i, d[seq].ratio)
                    elif self.running_options.get() == "label":
                        merge_dict[seq].add_all_label_mode(i, d[seq].counter, d[seq].peak_area,
                                                           d[seq].peak_intensity, d[seq].rt_seconds, d[seq].ions)

                else:
                    assert merge_dict[seq]
                    merge_dict[seq].add_sum_count(d[seq].counter)
                    if self.running_options.get() == "default" or self.running_options.get() == "lysine":
                        merge_dict[seq].add_ratio(i, d[seq].ratio)

                    if self.running_options.get() == "label":
                        merge_dict[seq].add_all_label_mode(i, d[seq].counter, d[seq].peak_area,
                                                           d[seq].peak_intensity, d[seq].rt_seconds, d[seq].ions)

                # merge_dict[seq].print_united()
        #################################
        #TO DO - this should be in Model!! not in controller!!
        unite_header = self.model.header_unite_create(n, self.running_options.get())
        self.model.xlsx_create(self.output_entry.get(), merge_dict, unite_header, self.running_options.get())
        ################################3
        self.error_message("Success", "Done!")
        # did not find non-blocking join in python :( if i dont use os exit the other thread just keep running forever
        #os._exit(0)
        exit(0)

    def show_frame(self, curr_frame):
        frame = self.frames[curr_frame]
        frame.tkraise()

    def error_message(self, error_name, message):
        messagebox.showinfo(error_name, message)

    def main(self):
        self.mainloop()

class StartPage(tk.Frame):
    def __init__(self, controller: Controller):
        tk.Frame.__init__(self, controller.container)
        self.controller = controller
        self.change_label = self._design()  # number of files chosen

    def _design(self):
        """ function return label to adjust
        """
        request_file = tk.Label(self, text="Select all pep.xml files", font=MAIN_FONT, fg="LightBlue4")
        request_file.grid(column=0, row=0, padx=20, pady=10)
        button_browse = ttk.Button(self, text=" Browse ", command=self.controller.choose_file)
        button_browse.grid(column=0, row=1)
        file_ok = tk.Label(self, text="No files selected yet", font=("David", 10), fg="dim gray")
        file_ok.grid(column=0, row=2)
        request_output_name = tk.Label(self, text="Enter output file name", font=MAIN_FONT, fg="LightBlue4")
        request_output_name.grid(column=0, row=3, pady=5)
        self.controller.output_entry = tk.Entry(self)
        self.controller.output_entry.grid(column=0, row=4)
        request_error_rate = tk.Label(self, text="Select desired error rate", font=MAIN_FONT, fg="LightBlue4")
        request_error_rate.grid(column=0, row=5, pady=5)
        error_0 = tk.Radiobutton(self, text="0.000", variable=self.controller.error_rate_button, value="0.0000")
        error_001 = tk.Radiobutton(self, text="0.001", variable=self.controller.error_rate_button, value="0.0010")
        error_005 = tk.Radiobutton(self, text="0.005", variable=self.controller.error_rate_button, value="0.0050")
        error_01 = tk.Radiobutton(self, text="0.010", variable=self.controller.error_rate_button, value="0.0100")
        error_0.grid(column=0, row=6)
        error_001.grid(column=0, row=7)
        error_005.grid(column=0, row=8)
        error_01.grid(column=0, row=9)
        request_mod = tk.Label(self, text="Select running mode", font=MAIN_FONT, fg="LightBlue4")
        request_mod.grid(column=0, row=10, pady=5)
        default_mod = tk.Radiobutton(self, text="Default                    ",
                                     variable=self.controller.running_options, value="default")
        label_free_mod = tk.Radiobutton(self, text="Label free                 ",
                                        variable=self.controller.running_options,
                                        value="label")
        lysin_uniform_mod = tk.Radiobutton(self, text="K & uniform n-term", variable=self.controller.running_options,
                                           value="lysine")
        default_mod.grid(column=0, row=11)
        label_free_mod.grid(column=0, row=12)
        lysin_uniform_mod.grid(column=0, row=13)

        #label_free = tk.Checkbutton(self, text="Label free mode", variable=self.controller.label_free,
                                #font=MAIN_FONT, fg="LightBlue4")
        #label_free.grid(column=0, row=10, pady=7)
        # TO DO!!!
        # i change it to "run" (used to be "next") button in the gui because im not using the second page
        # so it is actualy calling the "next" function and "next"  function call "run" function inside
        run_button = ttk.Button(self, text="Run", command=lambda: self.controller.next())
        run_button.grid(column=0, row=14, padx=20, pady=15)
        # lbl = ImageLabel(self)
        # lbl.grid(column=0)
        # lbl.load('phage.gif')
        return file_ok

    def browse(self):
        # Build a list of tuples for each file type the file dialog should display
        my_file_types = [('xml filles', '.xml'), ('all files', '.*')]
        files = askopenfilename(parent=self,
                               initialdir=os.getcwd(),
                               title="Please select a file:",
                               filetypes=my_file_types,
                                multiple=True)
        return list(files)

    def change_number_of_files_label(self, n):
        text = str(n) + " files have been selected"
        self.change_label = ttk.Label(self, text=text, font=("David", 10))
        self.change_label.grid(column=0, row=2)




# class SecondPage(tk.Frame):
#     def __init__(self, controller: Controller):
#         tk.Frame.__init__(self, controller.container)
#         self.controller = controller
#         self.lbl = ImageLabel(self)
#         #self.dict_check_buttons = dict()
#         for item in controller.files:
#             self.controller.swap_dict[item] = tk.IntVar()
#         request_messege = tk.Label(self, text="Select  the \"swapped ratio\" files", font=MAIN_FONT,
#                                    fg="LightBlue4")
#         default_lable = tk.Label(self, text= "All unselected files get ratio = light/heavy", font=("David", 10), fg="LightBlue4")
#         request_messege.grid(column=0, row=0, pady=5, padx=10, columnspan=2)
#         default_lable.grid(column=0, row=1, columnspan=2)
#         assert(controller.number_of_files == len(controller.files))
#         for i in range(controller.number_of_files):
#             check_button = tk.Checkbutton(self, text=controller.files[i],
#                                           variable=self.controller.swap_dict[controller.files[i]],
#                                           font=("Arial", 9))
#             check_button.grid(column=0, row=i + 2, columnspan=2, sticky="W")
#
#         curr_row = controller.number_of_files + 2
#         back_button = ttk.Button(self, text="Back", command=lambda: controller.back_to_start_page())
#         back_button.grid(row=curr_row, column=0)
#         run_button = ttk.Button(self, text="Run",
#                                 command=lambda: self._progress_bar(run_button, back_button, curr_row))
#         run_button.grid(row=curr_row, column=1, padx=10, pady=10)
#
#
#
#
#     def _progress_bar(self, run: ttk.Button, back: ttk.Button, curr_row):
#         run.configure(state=DISABLED)
#         back.configure(state=DISABLED)
#         self.lbl.grid(column=0, row=curr_row + 1, columnspan=2)
#         self.lbl.load('phage.gif')
#         precent = tk.Label(self, text="Running...please wait", font=("David", 10),
#                                  fg="LightBlue4")
#         precent.grid(column=0, row=curr_row+2, columnspan=2)
#         running_thread = threading.Thread(target=self.controller.run)
#         running_thread.start()





# def utility_func(i):
#     for k in range(i):
#         print(k)
#     exit(2)



# class ImageLabel(tk.Label):
#     """a label that displays images, and plays them if they are gifs"""
#     def load(self, im):
#         if isinstance(im, str):
#             im = Image.open(im)
#         self.loc = 0
#         self.frames = []
#
#         try:
#             for i in count(1):
#                 self.frames.append(ImageTk.PhotoImage(im.copy()))
#                 im.seek(i)
#         except EOFError:
#             pass
#
#         try:
#             self.delay = im.info['duration']
#         except:
#             self.delay = 100
#
#         if len(self.frames) == 1:
#             self.config(image=self.frames[0])
#         else:
#             self.next_frame()
#
#     def unload(self):
#         self.config(image=None)
#         self.frames = None
#
#     def next_frame(self):
#         if self.frames:
#             self.loc += 1
#             self.loc %= len(self.frames)
#             self.config(image=self.frames[self.loc])
#             self.after(self.delay, self.next_frame)





if __name__ == '__main__':
    window = Controller()
    window.main()

