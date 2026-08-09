// ==========================================================
// RS-FISH TASK-LIST BATCH (CSV DRIVEN — FINAL WORKING VERSION)
// Uses the SAME Bio-Formats import string as your original macro
// Supports 2- or 3-column CSVs
// Auto-detects channel, full path, base name
// Custom output directory
// ==========================================================
macro "RS-FISH Batch from CSV (FINAL — Original Import)" {

    // ================= USER CONFIG =================
    csvPath = "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/AP_file_2026/AP_files/rerun.csv";

    inputRoot  = "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/AP_file_2026/AP_files/";
    outBaseDir = "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/AP_file_2026/RSFISH_OUTPUT/";

    sigmaVal    = 1.8;
    anisotropy  = 0.81;
    maxError    = 1.5;
    inlierRatio = 0.10;
    supportVal  = 3;

    thr_CH1 = 0.0027863602;
    thr_CH2 = 0.007550505;
    thr_CH3 = 0.010730815;

    // -------------------------------------------------
    // Ensure output directory exists
    // -------------------------------------------------
    if (!File.exists(outBaseDir))
        File.makeDirectory(outBaseDir);

    // =================================================
    // READ CSV
    // =================================================
    csvText = File.openAsString(csvPath);
    if (csvText == "")
        exit("❌ Could not read CSV: " + csvPath);

    lines = split(csvText, "\n");

    print("Loaded task list: " + csvPath);
    print("Total tasks: " + (lines.length - 1));

    // =================================================
    // MAIN LOOP
    // =================================================
    for (i = 1; i < lines.length; i++) {

        line = trim(lines[i]);
        if (line == "") continue;

        cols = split(line, ",");

        baseName = "";
        imgPath  = "";
        chLabel  = "";

        // -------------------------------------------------
        // Detect channel (any column)
        // -------------------------------------------------
        for (c = 0; c < cols.length; c++) {
            token = trim(cols[c]);
            if (startsWith(token, "CH")) {
                chLabel = token;
                break;
            }
        }

        if (chLabel == "") {
            print("❌ Could not detect channel in line:");
            print(line);
            continue;
        }

        // -------------------------------------------------
        // Detect full path (any column containing / or \ )
        // (Your CSV has backslashes sometimes)
        // -------------------------------------------------
        for (c = 0; c < cols.length; c++) {
            token = trim(cols[c]);
            if (indexOf(token, "/") >= 0 || indexOf(token, "\\") >= 0) {
                imgPath = token;
                break;
            }
        }

        // -------------------------------------------------
        // If no full path, use base name (first column)
        // -------------------------------------------------
        if (imgPath == "") {
            baseName = trim(cols[0]);
            imgPath = baseName;
        }

        // -------------------------------------------------
        // Ensure .czi extension
        // -------------------------------------------------
        if (!endsWith(imgPath, ".czi") && !endsWith(imgPath, ".CZI"))
            imgPath = imgPath + ".czi";

        // -------------------------------------------------
        // Reconstruct full path if needed
        // -------------------------------------------------
        if (indexOf(imgPath, "/") < 0 && indexOf(imgPath, "\\") < 0) {
            partsTmp = split(imgPath, "_");
            geneTmp  = partsTmp[0];
            imgPath  = inputRoot + geneTmp + File.separator + imgPath;
        }

        // -------------------------------------------------
        // CHANNEL MAP (exact labels)
        // -------------------------------------------------
        if (chLabel == "CH1_633nm") {
            cIndex = 1; threshDoG = thr_CH1;
        } else if (chLabel == "CH2_555nm") {
            cIndex = 2; threshDoG = thr_CH2;
        } else if (chLabel == "CH3_488nm") {
            cIndex = 3; threshDoG = thr_CH3;
        } else {
            print("❌ Unknown channel: " + chLabel);
            continue;
        }

        print("\n--------------------------------------");
        print("PROCESSING:");
        print("Image  : " + imgPath);
        print("Channel: " + chLabel);

        if (!File.exists(imgPath)) {
            print("❌ IMAGE FILE NOT FOUND — SKIPPING");
            continue;
        }

        // -------------------------------------------------
        // Cleanup before processing
        // -------------------------------------------------
        run("Close All");
        winList = getList("window.titles");
        for (w = 0; w < winList.length; w++) {
            selectWindow(winList[w]);
            run("Close");
        }
        run("Collect Garbage");

        // -------------------------------------------------
        // Output directory (grouped by gene)
        // -------------------------------------------------
        fileName = File.getName(imgPath);
        dot = lastIndexOf(fileName, ".");
        base = substring(fileName, 0, dot);

        parts = split(fileName, "_");
        geneName = parts[0];

        outDir = outBaseDir + geneName + File.separator;
        if (!File.exists(outDir))
            File.makeDirectory(outDir);

        // =========================================================
        // ✅ BIO-FORMATS IMPORTER — EXACTLY LIKE YOUR ORIGINAL MACRO
        // =========================================================
        run("Bio-Formats Importer",
            "open=[" + imgPath + "] autoscale color_mode=Default " +
            "rois_import=[ROI manager] specify_range " +
            "view=Hyperstack stack_order=XYCZT series=1 " +
            "c_begin=" + cIndex + " c_end=" + cIndex +
            " z_begin=0 z_step=1");

        wImg = getWidth();
        hImg = getHeight();
        makeRectangle(0, 0, wImg, hImg);
        imgTitle = getTitle();

        // -------------------------------------------------
        // RUN RS-FISH (same as original)
        // -------------------------------------------------
        run("RS-FISH",
            "image=[" + imgTitle + "] mode=Advanced " +
            "anisotropy=" + anisotropy +
            " robust_fitting=RANSAC compute_min/max use_anisotropy " +
            "spot_intensity=[Linear Interpolation] " +
            "sigma=" + sigmaVal +
            " threshold=" + threshDoG +
            " support=" + supportVal +
            " min_inlier_ratio=" + inlierRatio +
            " max_error=" + maxError +
            " spot_intensity_threshold=0 background=[No background subtraction]");

        spots = nResults;
        print("Spots found: " + spots);

        // -------------------------------------------------
        // SAVE OUTPUTS
        // -------------------------------------------------
        outCSV = outDir + base + "_" + chLabel + "_RSFISH.csv";
        saveAs("Results", outCSV);

        fullLog = getInfo("log");
        outLog = outDir + base + "_" + chLabel + "_RSFISH_LOG.txt";
        File.saveString(fullLog, outLog);

        print("Saved CSV: " + outCSV);
        print("Saved LOG: " + outLog);

        // Cleanup after each task
        run("Close All");
        tableWins = getList("window.titles");
        for (t = 0; t < tableWins.length; t++) {
            if (tableWins[t] != "Log") {
                selectWindow(tableWins[t]);
                run("Close");
            }
        }
        run("Collect Garbage");
    }

    print("\n==== ALL TASKS COMPLETED SUCCESSFULLY ====");
}
