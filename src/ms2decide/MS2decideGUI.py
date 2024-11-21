import tkinter as tk
from tkinter import filedialog, messagebox
from pathlib import Path
from .K_estimation_GUI import K_estimation_GUI

# Your workflow function (K_estimation_GUI) is assumed to be imported and ready for use.

class MS2decideGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("K Estimation Workflow")
        
        # Define variables for the paths and values
        self.mgf_path = tk.StringVar()
        self.quan_path = tk.StringVar()
        self.job_description = tk.StringVar()
        self.isdb_path = tk.StringVar()
        self.ion_mode = tk.StringVar(value="Pos")  # Default to Pos
        self.tolerance = tk.DoubleVar(value=0.02)
        self.sirius_path = tk.StringVar()
        self.sirius_score_index = tk.StringVar(value="Exact")  # Default to Exact
        self.save_path = tk.StringVar()
        
        self.setup_ui()

    def setup_ui(self):
        # Authentication Section - Top Left
        auth_frame = tk.Frame(self)
        auth_frame.grid(row=0, column=0, padx=20, pady=20, sticky="nw")

        tk.Label(auth_frame, text="Authentication Section").grid(row=0, column=0, columnspan=2)

        tk.Label(auth_frame, text="Username:").grid(row=1, column=0)
        self.username_entry = tk.Entry(auth_frame)
        self.username_entry.grid(row=1, column=1)

        tk.Label(auth_frame, text="Password:").grid(row=2, column=0)
        self.password_entry = tk.Entry(auth_frame, show="*")
        self.password_entry.grid(row=2, column=1)

        tk.Label(auth_frame, text="Mail:").grid(row=3, column=0)
        self.mail_entry = tk.Entry(auth_frame)
        self.mail_entry.grid(row=3, column=1)

        # GNPS Section - Top Right
        gnps_frame = tk.Frame(self)
        gnps_frame.grid(row=0, column=1, padx=20, pady=20, sticky="ne")

        tk.Label(gnps_frame, text="GNPS Section").grid(row=0, column=0, columnspan=2)

        tk.Label(gnps_frame, text="MGF Path:").grid(row=1, column=0)
        tk.Entry(gnps_frame, textvariable=self.mgf_path).grid(row=1, column=1)
        tk.Button(gnps_frame, text="Browse", command=self.select_mgf_path).grid(row=1, column=2)

        tk.Label(gnps_frame, text="Quantitative CSV Path:").grid(row=2, column=0)
        tk.Entry(gnps_frame, textvariable=self.quan_path).grid(row=2, column=1)
        tk.Button(gnps_frame, text="Browse", command=self.select_quan_path).grid(row=2, column=2)

        tk.Label(gnps_frame, text="Job Description:").grid(row=3, column=0)
        tk.Entry(gnps_frame, textvariable=self.job_description).grid(row=3, column=1)

        # ISDB Section - Bottom Left
        isdb_frame = tk.Frame(self)
        isdb_frame.grid(row=1, column=0, padx=20, pady=20, sticky="sw")

        tk.Label(isdb_frame, text="ISDB Section").grid(row=0, column=0, columnspan=2)

        tk.Label(isdb_frame, text="ISDB Path:").grid(row=1, column=0)
        tk.Entry(isdb_frame, textvariable=self.isdb_path).grid(row=1, column=1)
        tk.Button(isdb_frame, text="Browse", command=self.select_isdb_path).grid(row=1, column=2)

        tk.Label(isdb_frame, text="Ion Mode:").grid(row=2, column=0)
        ion_mode_dropdown = tk.OptionMenu(isdb_frame, self.ion_mode, "Pos", "Neg")
        ion_mode_dropdown.grid(row=2, column=1)

        tk.Label(isdb_frame, text="Tolerance:").grid(row=3, column=0)
        tk.Entry(isdb_frame, textvariable=self.tolerance).grid(row=3, column=1)

        # Sirius Section - Bottom Right
        sirius_frame = tk.Frame(self)
        sirius_frame.grid(row=1, column=1, padx=20, pady=20, sticky="se")

        tk.Label(sirius_frame, text="Sirius Section").grid(row=0, column=0, columnspan=2)

        tk.Label(sirius_frame, text="Sirius Path:").grid(row=1, column=0)
        tk.Entry(sirius_frame, textvariable=self.sirius_path).grid(row=1, column=1)
        tk.Button(sirius_frame, text="Browse", command=self.select_sirius_path).grid(row=1, column=2)

        tk.Label(sirius_frame, text="Sirius Score Index:").grid(row=2, column=0)
        score_index_dropdown = tk.OptionMenu(sirius_frame, self.sirius_score_index, "Exact", "Approximate")
        score_index_dropdown.grid(row=2, column=1)

        # Save Path Button
        save_path_frame = tk.Frame(self)
        save_path_frame.grid(row=2, column=0, columnspan=2, padx=20, pady=20)

        tk.Button(save_path_frame, text="Select Save Path", command=self.select_save_path).grid(row=0, column=0)

        # Launch Workflow Button
        launch_button_frame = tk.Frame(self)
        launch_button_frame.grid(row=3, column=0, columnspan=2, padx=20, pady=20)

        tk.Button(launch_button_frame, text="Launch Workflow", command=self.launch_workflow).grid(row=0, column=0)

    # File Selection Methods
    def select_mgf_path(self):
        path = filedialog.askopenfilename(filetypes=[("MGF Files", "*.mgf")])
        if path:
            self.mgf_path.set(path)

    def select_quan_path(self):
        path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
        if path:
            self.quan_path.set(path)

    def select_isdb_path(self):
        path = filedialog.askopenfilename(filetypes=[("ISDB Files", "*.mgf")])
        if path:
            self.isdb_path.set(path)

    def select_sirius_path(self):
        path = filedialog.askopenfilename(filetypes=[("Sirius Files", "*.tsv")])
        if path:
            self.sirius_path.set(path)

    def select_save_path(self):
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.tsv")])
        if path:
            self.save_path.set(path)

    def launch_workflow(self):
        # Get the values from the UI and pass them to the K_estimation_GUI function
        username = self.username_entry.get()
        password = self.password_entry.get()
        mail = self.mail_entry.get()
        mgf_path = self.mgf_path.get()
        quan_path = self.quan_path.get()
        job_description = self.job_description.get()
        isdb_path = self.isdb_path.get()
        ion_mode = self.ion_mode.get()
        tolerance = self.tolerance.get()
        sirius_path = self.sirius_path.get()
        sirius_score_index = self.sirius_score_index.get()
        save_path = self.save_path.get()

        try:
            # Assuming K_estimation_GUI is imported correctly and available
            result = K_estimation_GUI(username, password, mail, mgf_path, quan_path, job_description,
                                      isdb_path, ion_mode, tolerance, sirius_path, sirius_score_index,save_path)
            messagebox.showinfo("Success", "Workflow completed successfully!")
            print(result)  # You can replace this with any output handling you need.
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {e}")

# Running the application
if __name__ == "__main__":
    app = K_estimation_GUI()
    app.mainloop()
