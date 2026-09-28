# SatelliteQ V26 — complete installation and user guide

This guide explains how to go from an unprocessed microscopy image to a Fiji segmentation result, a consolidated Python Excel report and interactive graphs. It is written for users who have basic familiarity with Fiji but little or no Python experience. **The installation steps are independent:** installing Python libraries does not install Fiji or its plugins.

**Start here:** install Fiji and Python, run **one single image in Fiji**, process the resulting Fiji file in Python, and inspect its Excel report and spatial viewer. Only after that should you batch or combine entire experiments.

## Contents

- [Part A — Installation and checking the computer](#part-a-installation-and-checking-the-computer)
- [Part B — Fiji step by step](#part-b-fiji-step-by-step)
- [Part C — Python step by step](#part-c-python-step-by-step)
- [Part D — What the measurements mean](#part-d-what-the-measurements-mean)
- [Part E — Experimental design and reproducibility](#part-e-experimental-design-and-reproducibility)
- [Part F — Troubleshooting](#part-f-troubleshooting)
- [Part G — What to preserve for a publication](#part-g-what-to-preserve-for-a-publication)

---

## Part A — Installation and checking the computer

### A1. What to download

1. [Fiji](https://fiji.sc/) for Windows, macOS or Linux.
2. `SatelliteQ_V26.ijm`, the Fiji macro.
3. Python 3.10 with working Tk/Tcl GUI support (the targeted development version; other Python versions should be tested before use).
4. `SatelliteQ_Python_V26.py`, the report generator.
5. `requirements.txt`, the five pip packages required by the Python script.
6. A modern web browser to display the generated HTML viewers.

Keep the two scripts at the repository filenames shown above. They do not have to be in the same folder as the microscopy images, but the Python script and `requirements.txt` should initially be together for easier installation.

### A2. Install and update Fiji

1. Download Fiji and extract or install it. Start Fiji.
2. In Fiji, choose **Help → Update…**. Install available updates and restart.
3. For **3D object analysis**, go to **Help → Update… → Manage update sites**.
4. Find and enable **3DObjectsCounterPlus**. If you do not see it, click **Add Unlisted Site** and supply the name `3DObjectsCounterPlus` and URL `https://sites.imagej.net/3DObjectsCounterPlus/`.
5. Apply the updates and restart Fiji. Confirm the menu command **Analyze → 3D Objects Counter+** is available. This is not the original **Analyze → 3D Objects Counter**. Documentation: https://imagej.net/plugins/3d-objects-counter-plus
6. If you plan to quantify cilia, also check **Skeletonize (2D/3D)** and **Analyze Skeleton**. Fiji's regular Auto Local Threshold, Bio-Formats, ROI Manager and Analyze Particles facilities are used in the appropriate branches.

If 3D Objects Counter+ still does not appear, open Fiji's updater again, verify that the site is enabled and that all updates have been applied, and restart Fiji once more.

### A3. Windows: install Python and the packages

Install a Python version with Tkinter, such as Python 3.10 from python.org. On Windows, enabling **Add Python to PATH** during installation can help, but you can also use the Python launcher (`py`).

Open the extracted SatelliteQ repository folder in Windows Explorer. Click its address bar, type `powershell`, and press Enter. Enter the following commands **one at a time**:

```powershell
py -3.10 --version
py -3.10 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pip check
.\.venv\Scripts\python.exe -m tkinter
```

A small Tk test window should appear after the last command. Close it. Use the same environment to run SatelliteQ:

```powershell
.\.venv\Scripts\python.exe SatelliteQ_Python_V26.py
```

**Why use `.venv`?** A virtual environment keeps SatelliteQ's packages separate from the rest of your computer. You can invoke its `python.exe` directly, so you do not need to change PowerShell execution policies to activate it.

If `py -3.10` is unavailable, use the command for your installed Python version (e.g., `py -3.11`) or install Python 3.10. This does not imply that all versions have been fully validated with V26.

### A4. macOS: install Python and the packages

Install a Python distribution with Tkinter support. Check the available version in Terminal:

```bash
python3 --version
```

Open Terminal in the SatelliteQ repository folder (`cd` to that directory), then run:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip check
.venv/bin/python -m tkinter
```

Close the Tk test window. Start SatelliteQ with:

```bash
.venv/bin/python SatelliteQ_Python_V26.py
```

If `-m tkinter` fails on macOS, use an installation of Python bundled with functioning Tcl/Tk. On some systems the default system Python has no working graphical Tk support.

### A5. Linux note

The Python virtual-environment commands are identical to macOS. Your distribution may supply Tkinter separately (often as `python3-tk`). Use a graphical desktop session for the Fiji dialogs and Python file picker.

### A6. Check installed package versions

Within the virtual environment, verify the five required packages:

```bash
python -c "import numpy,pandas,matplotlib,plotly,openpyxl; print('Five packages available')"
```

If you did not activate the environment, use `.venv/bin/python` (macOS/Linux) or `.\.venv\Scripts\python.exe` (Windows) in place of `python`.

Record the exact dependency versions actually used for an experiment:

```bash
python -m pip freeze > installed_versions.txt
```

The version ranges in `requirements.txt` are proposed installation constraints, not proof that every possible combination has passed end-to-end Fiji/Python testing.

### A7. Computer and storage planning

An ordinary laptop or desktop is sufficient to get started, and no GPU is required by the Python package list. A multi-core CPU, approximately 16 GB RAM and free storage for raw data, segmentation images and reports are practical starting points for typical confocal datasets, **not minimum verified specifications**. Memory consumption depends heavily on XY image dimensions, channel number, Z-stack depth, frame count and how many Fiji image duplicates are temporarily opened.

For first-time use, analyze one image instead of an entire time-lapse or an entire folder.

---

## Part B — Fiji step by step

### B1. Prepare and inspect the raw image

1. Open your fluorescence image or hyperstack in Fiji using **File → Open…** or the Bio-Formats importer for your microscopy format.
2. Confirm the image is the intended cell/field and contains a satellite fluorescence channel.
3. Check the **channel order**. For example, PCM1 may be channel 1, centrosome marker channel 2, cytoplasmic or membrane marker channel 3 and cilium marker channel 4. Your actual channel numbers may differ.
4. Inspect **Image → Properties…**. Verify pixel width and height, spatial unit, number of Z-slices and Z-step. Do not confuse an image Z-series with channels or time frames.
5. Keep the raw image unchanged. SatelliteQ duplicates and processes channels during analysis.

**Required input:** satellite fluorescence channel. **Optional:** centrosome, membrane/cytoplasmic cell boundary and cilium channels.

### B2. Start SatelliteQ

With the image already open, choose **Plugins → Macros → Run…**, select `SatelliteQ_V26.ijm`, read the welcome window and click OK. Choose one of the four analysis modes:

- **2D single image:** XY measurement on an image or a projected stack; area is measured in µm².
- **3D single stack:** object segmentation across a full Z-stack; apparent volume is measured in µm³.
- **2D time-lapse:** repeats the 2D measurement for a user-selected time range.
- **3D time-lapse:** repeats the 3D measurement for a user-selected time range.

Choose the analysis mode **before** editing segmentation settings; some options and unit conventions differ between modes.

### B3. Choose channels and projection (2D)

The 2D single-image interface has four major groups of settings. In **Step 1/4: Channels**, enter the actual **1-based Fiji channel number** for the satellite marker. For optional channels, enter their number if present or use `0` when the dialog says that skipping is allowed:

- Satellite marker channel (required).
- Centrosome marker channel (optional; for distance and PCI).
- Membrane/cytoplasmic boundary channel (optional; for automatic cell boundaries).
- Primary cilium marker channel (optional; for cilium analysis).

If your input is a Z-stack, choose **Max Intensity**, **Sum Slices** or **Average Intensity** projection. The default projected analysis commonly uses Max Intensity. **Document the chosen projection and use the same method for directly compared conditions.** A projected image loses axial information: satellite and centrosome distances are XY distances.

### B4. Choose satellite settings

In **Step 2/4: Satellites**, you can set:

- **Rolling-ball radius:** controls subtraction of slowly varying background, in pixels. Choose a radius using representative images, not just the brightest field.
- **Threshold mode:** `Auto`, `Manual value` or `Hysteresis-like` (2D).
- **2D automatic methods:** `Otsu`, `Triangle`, `Moments`, `Default`, `Li` and `Yen`.
- **Manual threshold:** a user-entered fluorescence cutoff when `Manual value` is selected.
- **Hysteresis-like high-threshold multiplier:** used only when that mode is selected.
- **Minimum/maximum satellite size:** reject regions that are too small or too large for your experimental question; check the units in the mode-specific dialog and output settings.
- **Binary watershed:** optional separation of adjacent 2D objects.
- **Large-cluster splitting:** optional `Local intensity peaks` branch for large merged 2D regions, with settings for minimum cluster size, peak-centered spot radius and prominence.
- **Review/edit detected satellites:** enable when you want to inspect or correct the results before quantification.

**How to choose the settings:** On a representative image, look for a threshold that detects the visible satellite signal without turning cytoplasmic background into many small false-positive objects. Compare the segmentation overlay directly with the raw fluorescence. Check the effect of minimum object size on dim satellites and of watershed/large-cluster splitting on crowded pericentrosomal regions. Once chosen, apply identical settings to the compared control and experimental images from the same acquisition dataset.

Do not assume that using the identical numeric threshold on very different staining or microscope acquisitions will create equivalent detection sensitivity; document any dataset-specific optimization.

### B5. Centrosome detection and review

If a centrosome channel was provided, the 2D workflow offers **automatic centrosome detection**, **manual centrosome ROIs during the macro** or **existing centrosome ROIs**. For automatic detection you can adjust background subtraction, Gaussian smoothing, bright-signal thresholding, object-area limits and the maximum number of centrosomes retained. Review the resulting markers against the centrosome channel and satellite channel when possible.

The 3D workflow uses the centrosome channel to identify XY positions on a projection and estimates each centrosome's Z-position from the centrosome-channel Z-stack. Because that Z estimate can be sensitive to ROI placement, noise and local intensity, **verify representative centrosome XYZ positions against the original stack before relying on 3D distances or PCI**.

When multiple centrosomes are associated with the same cell, SatelliteQ offers policies such as assigning satellites to the nearest centrosome, using a midpoint reference or excluding the cell from centrosome-dependent analyses. Record the policy in your methods; the choices can change the resulting distance and PCI.

### B6. Cell assignment: choose the biological population you want

This step is particularly important: **a different cell-selection region changes which satellites contribute to the per-cell statistics**.

In 2D single-image mode, the available strategies include:

- **No cell assignment:** useful for image-level measurements; not suitable for a per-cell biological comparison without additional assignments.
- **Membrane marker:** derive cell ROIs from a membrane/boundary channel.
- **Cytoplasmic signal boundary:** derive cell outlines from a suitable cytoplasmic fluorescence channel.
- **Centrosome radius:** use defined centrosome-centered circular regions; not identical to whole-cell analysis.
- **Manual cell ROIs during macro:** draw and review cell boundaries yourself.
- **Existing cell ROIs in ROI Manager:** reuse previously prepared ROIs.

The 3D single-stack workflow includes options such as **centrosome-centered spheres**, projected 2D cell ROIs derived from membrane/cytoplasmic channels and true 3D cell-region options. The filtering units can differ between projected ROI areas and true 3D cell volumes; read the specific dialog labels carefully. Choose how to handle satellites outside a selected cell region (ignore them or export them as unassigned, when offered).

If you use automatic boundary detection, inspect the mask and cell labels. Local-threshold options include Phansalkar, Sauvola, Bernsen, Niblack, Mean and Median for the cell marker. Additional cleanup, preprocessing and watershed options are available for difficult markers such as cytoskeletal staining.

**Avoid changing the assignment mode between control and treatment** in a direct comparison. A centrosome-centered sphere may include a different proportion of peripheral satellites than a full cellular ROI.

### B7. Optional primary cilium channel

When a cilium marker channel is provided and cilium detection is enabled, set its background subtraction and smoothing, threshold method and minimum/maximum accepted path length. The default minimum accepted length in V26 is **1 µm**, and users can change it.

Cilium length is derived from processed **2D projected masks** and skeleton path measurements, even when the satellite input is 3D. Review cilium detections and check that nonciliary structures are not included. The optional cilium readouts can be associated with the same cell IDs as the satellite measurements when valid cell assignment is available.

### B8. Special points for 3D satellite analysis

Install and check **3D Objects Counter+** *before* starting a 3D job. Use an image with real Z metadata or enter the correct calibration when asked. The 3D workflow provides **Otsu** as its automatic satellite threshold or permits a manual threshold value; do not expect the entire list of 2D automatic satellite thresholds to be available in 3D. Review object-size limits and the option to exclude stack-edge objects when offered.

After a small trial analysis, audit several exported satellite coordinates against the source stack, including a few satellites near and far from the centrosome. Also verify at least one centrosome's estimated Z-slice. Do not treat an unexpectedly large 3D distance as a biological phenotype until coordinate units, Z-index convention and centrosome position have been checked.

### B9. Time-lapse

For 2D or 3D time-lapse input, choose the frame range and interval in the mode-specific dialogs, then check segmentation and assignment on representative early, middle and late frames before running the entire series. Save the time interval and relevant acquisition metadata.

Each frame is quantified separately. You can subsequently plot satellite number, total/mean size, fluorescence intensity and other available per-frame measurements. **V26 does not track individual satellite identities across frames** and cannot by itself prove that one specific satellite fused with another or moved along a particular path.

### B10. Export Fiji results and review them

The generated filenames vary by analysis mode and run. Keep the **original names** because the Python program identifies recognized patterns. For example, 2D results often end in `_SatelliteQuantify_Fiji_Output.csv`; 3D outputs often follow `ImageName_3D_<timestamp>_results.csv`; time-lapse outputs can include `timelapse_satellite_rows_combined.csv`.

Preserve available settings rows, calibration information, QC images, cell-boundary coordinate files, cilium companion tables and the original raw microscopy data. Confirm that one of the exported results files contains legitimate satellite rows before opening Python.

---

## Part C — Python step by step

### C1. Start the program

On Windows, from the repository folder in PowerShell:

```powershell
.\.venv\Scripts\python.exe SatelliteQ_Python_V26.py
```

On macOS/Linux, from the repository folder in Terminal:

```bash
.venv/bin/python SatelliteQ_Python_V26.py
```

The Tkinter dialog asks whether you wish to choose a **folder containing one or more SatelliteQ outputs** or **one CSV/Excel file**. You can bypass that dialog by giving the path directly:

```bash
python SatelliteQ_Python_V26.py "/path/to/the/Fiji_output.csv"
python SatelliteQ_Python_V26.py "/path/to/experiment"
```

If `python` launches the wrong interpreter, replace it with your `.venv` Python as above.

### C2. Start with one Fiji result file

1. Select the original result CSV/TSV/XLSX exported by Fiji (not the microscopy TIFF or a Python-generated report).
2. Wait for the program to read the satellite and related companion tables, standardize field names, calculate available summaries and generate plots/viewers.
3. Check the new `<your input name>_Python_Report/` folder beside that input file.
4. Open the generated quantification `.xlsx` report and inspect the `README`, `Input_files`, `Summary` and satellite/cell worksheets.
5. Open the Graph Explorer and Spatial Viewer HTML outputs. Compare a representative cell's count, intensity and position with the corresponding Fiji segmentation.

**If no usable satellite rows are detected**, first open the Fiji output itself and make sure it contains satellite objects and retained its original Fiji filename format. Python cannot compute cell summaries from images alone.

### C3. Combine conditions and biological replicates

Put the Fiji result files for a multi-condition experiment into folders that clearly identify the biological groups. Recommended layout:

```text
Experiment_01/
  Control/
    Replicate_1/
      <Fiji-generated output>.csv
    Replicate_2/
      <Fiji-generated output>.csv
  Nocodazole/
    Replicate_1/
      <Fiji-generated output>.csv
    Replicate_2/
      <Fiji-generated output>.csv
```

Choose `Experiment_01` when prompted to select a folder. The program searches recursively for recognized Fiji files and, if it finds several, lets you choose whether to **combine**, **process separately** or **do both**.

- **Combine** creates an experiment-wide `SatelliteQ_Combined_Report/` with condition and replicate labels.
- **Separate** creates individual `<input name>_Python_Report/` folders.
- **Both** produces both forms of output.

Keep generated report folders separate from source inputs and avoid renaming original Fiji files. The program is designed to ignore common previous report names, but an intentionally clean input folder makes auditing simpler.

### C4. Optional experiment metadata file

If your folder names do not reliably encode the experimental groups, add a comma-separated file named **`SatelliteQ_experiment_metadata.csv`** to the experiment folder, using this header and your **actual Fiji result filenames**:

```csv
Input_file,Condition,Replicate
Image_01_SatelliteQuantify_Fiji_Output.csv,Control,Replicate_1
Image_02_SatelliteQuantify_Fiji_Output.csv,Control,Replicate_2
Image_03_SatelliteQuantify_Fiji_Output.csv,Nocodazole,Replicate_1
Image_04_SatelliteQuantify_Fiji_Output.csv,Nocodazole,Replicate_2
```

The filenames above are **illustrative**: replace them with the actual names exported by Fiji. V26 searches for this metadata file (and recognized alternatives such as `experiment_metadata.csv` or `conditions.csv`) and uses matching rows to annotate Condition and Replicate. Keep metadata as a plain CSV file and confirm it matches the input filenames exactly.

### C5. What Python produces

For every successful analysis, V26 creates an Excel report that includes the available worksheets, which can include:

- `README`: report version and explanation.
- `Input_files`: input-file and path index.
- `Summary`: overall analysis summary.
- `Satellite_rows_clean` or `Combined_satellite_rows`: standardized object-level measurements.
- `Image_summary` or `Combined_image_summary`: measurements summarized by image.
- `Per_cell_summary` or `Combined_cell_summary`: one row per valid cell (and frame when time-lapse is present).
- `Spatial_metrics` or `Combined_spatial_metrics`: available spatial-organization metrics.
- `Frame_summary` or `Combined_frame_summary`: measurements summarized by time-lapse frame.
- `Cilium_rows`, `Cilium_cell_summary` and `Cilium_skeleton_points`, when present.
- `Units`: units and meaning of standardized columns.
- `Graph_and_viewer_files`: an index of files generated by the report.

The script can also write a separate explanatory Excel workbook. Not every output sheet appears for every experiment: the information required for a given analysis must be present in the Fiji input.

### C6. Graph Explorer and Experiment Viewer

V26 creates self-contained Plotly HTML viewers. The **General Graph Explorer** allows you to explore available quantitative variables and their relationships. The **Experiment Viewer** focuses on cell-level comparisons across experimental conditions and replicates, with display controls for the supported parameters. A **tabbed dashboard** can show the two viewers together.

Start with one biological variable (for example satellite number or mean centrosome distance), check that the labels and replicates are correct, then examine the other parameters. Use any normalization options deliberately and report the reference group in figure legends. Plots produced for exploration are not a substitute for selecting an appropriate statistical model that respects your experimental replicate structure.

### C7. Spatial Viewer

The Spatial Viewer displays satellites at their measured positions and can show centrosomes, cell boundaries and optional cilia when source data exist. Filter to a specific cell or frame, toggle layers and examine coloring/marker-size options to identify outliers or segmentation issues.

When Fiji provided exact 2D cell-boundary ROI coordinates, those outlines are displayed. Otherwise an **estimated boundary** may be drawn from the cell's available area/volume information. An estimated shape is a visualization aid, **not a replacement for the actual segmented cell boundary**.

Use the Spatial Viewer as a cross-check: a cell showing unusually high mean centrosome distance or PCI should be inspected against the raw microscopy image and Fiji ROI assignment.

### C8. File organization and practical workflow

For a publication, preserve separate folders for **raw images**, **original Fiji exports**, **Python reports** and **analysis-setting records**. Avoid editing CSV columns by hand to force analysis: change or correct segmentation in Fiji, then rerun Python. Keep a record of the V26 source-code revision and the exact Python package versions (`python -m pip freeze`).

---

## Part D — What the measurements mean

**Satellite number** is the count of accepted segmented objects. It depends on the acquisition quality, chosen threshold, object-size filters, adjacent-object separation and cell-selection region. A segmentation object may contain more than one tightly overlapping biological satellite.

**Area and volume** measure an object's fluorescent footprint in 2D and apparent fluorescence volume in 3D. Do not compare 2D area (µm²) directly with 3D volume (µm³). With diffraction-limited microscopy, these values are not electron-microscopy measurements of individual ultrastructural granule size.

**Integrated intensity** is the sum of background-corrected signal within an object. The 2D code uses Fiji `RawIntDen`. Intensity depends on acquisition and staining; comparisons should use matched settings and transparently described normalization where applicable.

**Distance to the centrosome** is the straight-line distance from satellite centroid to its assigned centrosome, using XY for 2D and XYZ for 3D, when valid coordinates exist. Distances depend on which satellites were included in cell assignment. If a centrosome is missing or its Z value is incorrect, the derived distance can be missing or misleading.

**Pericentrosomal clustering index (PCI)** is calculated from object-level integrated fluorescence: the total satellite intensity at centroid distance **0–3 µm** from the centrosome, divided by the total intensity at centroid distance **>3–12 µm**. The two radii are **fixed**, not adjustable, in V26. A zero outer-region intensity makes PCI **undefined**, not zero. In 2D the grouping uses XY radius; in 3D it uses XYZ radius.

**Spatial entropy** uses intensity weights assigned to XY grid locations of satellite centroids, at **5 × 5**, **10 × 10** and **20 × 20** resolutions. Each grid's Shannon entropy is normalized by `ln(min(K, N))`, and the reported value is the mean of the valid scales. Low values correspond to spatially concentrated XY signal; high values indicate broader or more even distribution across sampled XY bins. **The current 3D image workflow does not compute volumetric 3D entropy**: its segmentation is 3D, but its entropy binning remains XY.

**Shape** is measured as 2D circularity or 3D sphericity when the required object measurements are available. Tiny near-diffraction-limited satellites can yield unreliable shape values; interpretation is more straightforward for enlarged granules.

**Cilium analysis** uses projected/skeleton-based length measurements and a configurable ciliation cutoff (default minimum 1 µm). Even when the input stack is 3D, the cilium length measurement remains projected rather than the actual 3D axoneme length.

**Time-lapse summaries** describe how the satellite population changes within a cell over time. They do not follow object identity between frames.

---

## Part E — Experimental design and reproducibility

1. **Design the acquisition first.** Keep laser power, detector gain, pixel size, Z-step, optical zoom and fixation/staining protocols as comparable as possible between groups intended for direct comparison.
2. **Optimize on representative images.** Tune threshold and size criteria by visually reviewing representative control and treated cells with a range of signal intensities, not just the highest-quality image.
3. **Match settings within experiments.** Fix the analysis mode, projection method, thresholding, rolling-ball, size filters and cell-assignment policy for direct biological comparisons. Different experiments with different imaging/staining quality may require separate, documented optimization.
4. **Check cell assignment.** Whole-cell boundaries and centrosome-centered regions can include substantially different satellite populations and produce different per-cell counts and mean distances.
5. **Check 3D calibration.** Audit Z-step, coordinate conventions, centrosome axial position and exported measurements on a representative Z-stack. When possible, independently recalculate a few satellite-to-centrosome distances from verified XYZ coordinates.
6. **Check the source rows for every condition.** Compare satellite counts, intensity distributions, fraction of excluded objects, number of valid cells, and whether unusually high distances arise from true peripheral satellites or measurement/assignment errors.
7. **Analyze biological replicates appropriately.** Keep per-cell or per-granule rows linked to their original image, frame, condition and biological replicate. Do not present many satellite objects or cells as independent biological experiments.
8. **Keep a settings record.** Save Fiji settings/QC, Python reports, code commit and dependency snapshot. Show parameter selection and any major re-optimization transparently in the associated manuscript.

---

## Part F — Troubleshooting

**Fiji says no image is open.** Open the image/hyperstack before running the macro; do not launch the macro from an empty Fiji session.

**Channel numbers do not match.** Check the channel slider and actual channel order after Bio-Formats import. Reopen a single field to verify the image metadata; don't guess based on display colors.

**3D analysis cannot find 3D Objects Counter+.** Install the correct **3DObjectsCounterPlus** update site and restart Fiji; verify the command in the Analyze menu. The original 3D Objects Counter and the StarDist variant are different plugins.

**Segmentation detects too many tiny dots.** Review background subtraction and threshold selection, then review the minimum object size and the raw fluorescence. If results look different between experiments with very different image quality, do not assume identical numeric threshold values ensure equivalent detection.

**Segmentation merges nearby satellites.** Review the 2D watershed/large-cluster splitting options and the corresponding overlays. In 3D, inspect Z-slices directly; no algorithm can resolve individual satellites below the effective optical resolution of the input image.

**Cell counts/distances differ between figures.** Check whole-cell versus centrosome-centered regions, projected versus 3D cell-boundary handling, image quality and acquisition scale. Verify that control and treatment were analyzed with matched rules.

**3D distance or PCI is unexpectedly large.** First verify the physical Z-step; inspect a few satellite and centrosome XYZ coordinates in the raw stack; audit whether the centrosome Z estimate is plausible; then check which satellites were included in the analyzed cell. Do not automatically attribute outliers to biological dispersal.

**Python command not found / missing package.** Use the Python executable in `.venv`. Check `python -m pip check` and test imports with the same interpreter. Installing packages into a different Python installation will not fix the active environment.

**Tk dialog does not appear.** Test `python -m tkinter`; install a Tk-enabled Python distribution if needed. Alternatively, launch V26 with the original Fiji result file or parent folder as a command-line argument.

**Python finds no compatible files.** Point it at the original Fiji result file, or at a folder with recognized `*_SatelliteQuantify_Fiji_Output.csv`, 3D `*_3D_*_results.csv` or time-lapse `*timelapse_satellite_rows_combined.csv` names. Avoid manually renaming source data or supplying raw image files.

**Incorrect condition or replicate labels.** Put the files under `Experiment/Condition/Replicate/` or supply a correctly named `SatelliteQ_experiment_metadata.csv`. Check that its `Input_file` entries match the real filenames. Inspect the `Input_files` and summary sheets after the Python run.

**Spatial metric is blank.** Missing or invalid cell assignment, centrosome position, satellite distance or integrated intensity can prevent some calculations. For PCI, an empty outer 3–12 µm shell also makes the index undefined. Never replace missing metric values with zero without a documented scientific reason.

**HTML viewer looks different from the Fiji mask.** Check whether exact cell-boundary ROI coordinates were exported; the viewer may display an estimated boundary if no exact outline is available. The original Fiji image/mask remains the reference for segmentation.

---

## Part G — What to preserve for a publication

Archive the raw images, acquisition metadata, final Fiji macro version, original Fiji CSV/TSV and settings/QC images, final Python script revision, `installed_versions.txt`, input/metadata folder structure, unmodified Python Excel reports and viewer outputs. Add the agreed code license and citation metadata to the public repository when the release is finalized. A reproducible biological comparison should make clear the number of independent biological replicates, analyzed cells, image dimensions, chosen thresholds, cell-selection policy and any exclusions.
