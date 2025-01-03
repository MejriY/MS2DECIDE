## Usage Guide: `K_estimation`

The `K_estimation` function is a core feature of `ms2decide`, designed to aggregate annotations from multiannotation tools as GNPS, ISDB-LOTUS, and Sirius by leveraging decision theory. This approach models expert preference elicitation to derive \( K \)-values for features prioritization, providing a unified and informed annotation estimate.

---

### **Overview**
The `K_estimation` function:
- Integrates annotation data from multiple sources: GNPS, ISDB-LOTUS, and Sirius.
- Processes user-provided input files such as quantitative data and MGF files.
- Generates a filtered dataframe based on the estimated \( K \)-values.
- Exports the results into a `.tsv` file for further analysis.

---

### **How to Use `K_estimation`**

1. **Input credentials**
   - The function will prompt you to provide your GNPS username, password, and email. This is necessary for authenticating on the GNPS platform.

2. **Provide Input Files**
   - Specify the paths to the following required files:
     - **Quantitative Table File**: A `.csv` file containing the quantitative data for your analysis. Ensure the file format is correct to avoid errors.
     - **MGF File**: A file containing tandem mass spectrometry data in MGF format.
 
   
   - Example paths:
     ```plaintext
     SELECT THE PATH FOR YOUR QUANTITATIVE FILE. This path needs to terminate with a *.csv at the end 
     : /path/to/quantitative_file.csv

     SELECT THE PATH FOR YOUR MGF FILE.  This path needs to terminate with a *.mgf at the end 
     : /path/to/mgf_file.mgf
     ```

3. **Set Up GNPS Workflow**
   - **FBMN Title:** Provide a unique and descriptive title for your FBMN job. **Note**: A folder with this title will be created to upload the quantitative data and MGF files. Ensure the title you choose does not match the name of an existing folder.
   - Choose the type of GNPS library search:
     - `strict`: Uses a typical mass difference tolerance of 0.02 Da. For more precise information, we recommend visiting the [GNPS Documentation on Library Search](https://ccms-ucsd.github.io/GNPSDocumentation/librarysearch/). **Note**: This workflow uses a library score threshold of `0.001`.
     - `iterative`: for iterative weighted analog search (can take up to three hours).
At this level, 27 FBMN jobs will be launched on your GNPS account. In the case of `strict`, only one job will be launched.


image::https://github.com/MejriY/MS2DECIDE_pic/raw/main/image/gnps_iterative.png[]

4. **ISDB-LOTUS Annotation**
   - The ISDB-LOTUS annotation is performed using the function `isdb_res = get_cfm_annotation(mgf, ISDBtol)`. During the process, the user will be prompted to provide:
     - **Ionization Mode**: Specify the ionization mode for annotation (`POS` for positive, `NEG` for negative).
     - **Mass Tolerance**: Provide a mass tolerance value less than `0.5` (default: `0.02`). **Note**: This value is comprised between 0 and 0.5.
   - This function calculates annotations by matching mass spectrometry data against ISDB-LOTUS spectral data.

5. **Sirius Annotation**
   - Provide the path to the Sirius6 annotation file (`structure_identifications.tsv`).
   - Select the confidence score type:
     - `exact`
     - `approximate`

6. **Compile Annotations**
   - Annotations from GNPS, Sirius, and ISDB-LOTUS are compiled into a unified dataframe.
   - The dataframe is filtered and sorted by \( K \)-values.

7. **Export Results**
   - Specify the path to save the output `.tsv` file:
     ```plaintext
     SELECT THE SAVE PATH FOR THE .TSV FILE OF MS2DECIDE OUTPUT. 
     #This path needs to terminate with a file_name.tsv where `file_name` is the desired name specified by the user.
     ```

8. **Optional: Retrieve Empty Annotations in the case of iterative weighted GNPS analog search**
   - If requested (`yes`), the function generates a report of empty annotations and saves it as `empty.tsv`.

---

### **Return Value**
The function returns a (`tsv file`)containing the **processed** and **ranked** results.

---

### **Example Workflow**
Below is an example of how to call the `K_estimation` function in Python and interact with its prompts:

```python
from ms2decide.K_estimation import K_estimation

# Run the function to estimate K values
df = K_estimation()

# Output dataframe and file saved at specified path
print(df.head())
```

---

### **Notes**
- Ensure that your input files follow the required format (`.csv` for quantitative data, `.mgf` for mass spectrometry data).
- Ensure internet connectivity for accessing GNPS and related services.
- The save path for the `.tsv` file must have a valid `.tsv` extension and specify the desired file name (e.g., `/path/to/output_file.tsv`).

By following these steps, you can effectively use the `K_estimation` function to process and analyze your compound data.
