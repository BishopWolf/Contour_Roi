
import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.tool.PlugInTool;
import ij.plugin.filter.ThresholdToSelection;
import ij.plugin.frame.Recorder;
import ij.measure.Calibration;
import java.awt.*;
import java.awt.event.*;

/**
 * An ImageJ interactive magic wand with selectable tolerance, variable hue or 
 * grayscale, preference for RGB, gradient detection for grayscale, 4-connected,
 * 8-connected or disconnected operation and preview.
 *
 *
 * Parameters ==========
 *
 * Value Tolerance: The selection is expanded to all image points as long as the
 * difference between the pixel value of the point clicked and the image point
 * is less than the Value Tolerance. Ignored for thresholded images; there the
 * thresholded (red) area is selected instead. For 16-bit and float images, note
 * that rounding errors occur, thus the red area displayed is only an
 * approximation to the exact thresholded area. In contrast to the ImageJ wand,
 * you must click INTO the thresholded area; clicking outside will create a
 * beep. For color images, see the color slider below. The selection can be
 * further restricted by the Gradient Tolerance. For calibrated grayscale
 * images, the Gray Value Tolerance is in calibrated units (e.g., 'lux').
 *
 * Color: For RGB images, this slider determines the relative weight of the gray
 * value and the color hue/saturation for the evaluation of the tolerance. With
 * a value of -100, only the gray value (brightness) is taken into account (with
 * weights of the color channels depending on the
 * Edit>Options>Conversions>Weighted RGB Conversions checkbox). With a value of
 * 0, the vector distance between the (unweighted) RGB vectors is used, and with
 * a value of 100, only the difference of hue/saturation is taken into account,
 * independent of the gray value. Use a large value (e.g. 95) if you want to
 * select one color independent of brightness.
 *
 * Gradient Tolerance: Irrespective of the Value Tolerance and color settings,
 * the selection is not expanded if the gray level gradient is larger than the
 * Gradient Tolerance. Gradient tolerance is also active if a threshold is used
 * instead of the Value Tolerance; in this case visual inspection is easier with
 * the Threshold display in 'Over/Under' mode. Gradient detection is less
 * sensitive to noise than in 4-connected mode than in 8-connected mode. For
 * calibrated images, the gradient tolerance is in calibrated units, e.g.,
 * 'lux/mm'. For RGB images, gradient tolerance works on the gray value only.
 *
 * Eyedropper Color: Uses the ccurrent foreground color (color of the eyedropper
 * tool) as reference color (or pixel value), instead of the pixel color (value)
 * where the user has clicked. In 4-connected and 8-connected mode (see
 * 'connectedness', below) the starting pixel is selected irrespective of its
 * color (value). In non-contiguous mode, an empty selection may result if no
 * pixel matches the eydropper color within the tolerance.
 *
 * Connectedness: Can be "Eight-Connected", "Four_Connected" or
 * "Non-Contiguous": In Eight-Connected mode, the selection proceeds to
 * neighboring pixels touching each other at a side or a corner. Four-Connected
 * means that the selection proceeds only to neighboring pixels sharing a side
 * with the previous pixel. Non-Contiguous mode selects all pixels of the image
 * within the tolerance irrespective of connectedness. In non-contiguous mode,
 * the 'Gradient Tolerance' and 'Include Holes' options are ignored.
 *
 * Include Holes: When checking this option, the selection created by this wand
 * operation will not have interior holes. Nevertheless, holes may arise from
 * combining/subtracting selections with the SHIFT or ALT keys.
 *
 * Usage ===== This file should be in ImageJ/plugins/tools directory, then use
 * Compile&Run, then restart ImageJ (Fiji users currently need the .class files
 * or a .jar file). [For creating a .jar archive from .class files: jar -cvfM
 * Versatile_Wand_Tool.jar Versatile_Wand_Tool.*] Left-click the tool icon for
 * selecting the tool. Right-click or double-click the tool icon for the options
 * menu where the tolerance and other parameters can be selected. Note: For
 * large images, it is recommended to use ImageJ 1.50e or later.
 *
 *
 *
 * Version 1.21, Wayne Rasband, Michael Schmid 2015-Oct-04: Converted to
 * PlugInTool, fixed bug of Color sensitivity Version 1.15, Michael Schmid,
 * 2014-Jul-30: Bugfix, NullpointerException without image after 'do Wand on OK'
 * Version 1.14, Michael Schmid, 2012-Dec-09: Again Bugfix: RGB tolerance was
 * inconsistent Version 1.13, Michael Schmid, 2012-Dec-07: Bugfix RGB tolerance
 * was incorrect (incomplete) Version 1.12, Michael Schmid, 2010-Sep-10: Can use
 * eyedropper (foreground) color instead of the color (value) of the pixel
 * Version 1.11, Michael Schmid, 2009-Jul-10: Color, non-contiguous added.
 * NullPointerException fixed.
 *
 */
public class Contour_Roi_Tool extends PlugInTool implements DialogListener {

    private final static String CLASS_NAME = "Contour_Roi_Tool";
    private final static String[] CONNECTED_STRINGS = new String[]{"8-connected", "4-connected", "non-contiguous"};
    private final static int EIGHT_CONNECTED = 0, FOUR_CONNECTED = 1, NON_CONTIGUOUS = 2;
    private final static int UNKNOWN = 0, OUTSIDE = 1, INSIDE = -1;  //mask pixel values
    // dialog parameters
    private double valueTolerance = 0, valueTolerance0 = 0;
    private double gradientTolerance = 0;
    private double colorSensitivity = 0;
    private int connection = EIGHT_CONNECTED;
    private boolean includeHoles = true;
    private boolean useEyedropper = false;
    // remember previous position for preview
    private int xStart = -1, yStart = -1, yshift = 0;
    // preview thread
    private Thread previewThread;
    private boolean previewing;

    public void showOptionsDialog() {
        ImagePlus imp = WindowManager.getCurrentImage();
        int impType = 0;
        Roi originalRoi = null;
        if (imp != null && xStart < 0 || yStart < 0 || xStart >= imp.getWidth() || yStart > imp.getHeight()) {
            imp = null;                 // no starting position for preview, thus no preview image needed
        }
        if (imp != null) {
            if (imp.lockSilently()) // try locking the image, to avoid previewing with a locked image
            {
                imp.unlock();
            } else {
                IJ.showStatus("No preview (" + imp.getShortTitle() + " locked)");
                imp = null;
            }
        }
        if (imp != null) {              // we have an image for preview
            impType = imp.getType();
            originalRoi = imp.getRoi();
            imp.saveRoi();
        }
        String grayUnit = "", gradientUnit = "";
        if (imp != null && imp.getCalibration() != null) {
            Calibration cal = imp.getCalibration();
            grayUnit = "(" + cal.getValueUnit();
            gradientUnit = grayUnit + "/" + cal.getUnit() + ")";
            grayUnit = grayUnit + ")";
        }

        double saveValueTolerance = valueTolerance;
        double saveGradientTolerance = gradientTolerance;
        double saveColorSensitivity = colorSensitivity;
        int saveConnection = connection;
        boolean saveIncludeHoles = includeHoles;
        boolean saveUseEyedropper = useEyedropper;

        GenericDialog gd = new GenericDialog("Contour Roi Options");
        if (imp == null || impType == ImagePlus.COLOR_RGB
                || (impType == ImagePlus.GRAY8 && imp.getProcessor().getCalibrationTable() == null)) {
            gd.addSlider("Value Tolerance", 0, 255, valueTolerance);
        } else {
            gd.addNumericField("Value Tolerance", valueTolerance, 2, 6, grayUnit);
        }
        gd.addSlider("Color: Gray only - Std - Hue only", -100, +100, 100 * colorSensitivity);
        gd.addNumericField("Gradient Tolerance (0=off)", gradientTolerance, 2, 6, gradientUnit);
        gd.addCheckbox("Eyedropper Color (instead of pixel color)", useEyedropper);
        gd.addChoice("Connectedness", CONNECTED_STRINGS, CONNECTED_STRINGS[connection]);
        gd.addCheckbox("Include Holes", includeHoles);
        gd.addHelp("http://imagejdocu.tudor.lu/doku.php?id=plugin:segmentation:versatile_wand:start&dummy=.html");
        gd.addDialogListener(this);     //the DialogItemChanged method will be called on user input
        if (imp != null) {
            startPreview(imp, xStart, yStart);
        }
        gd.showDialog();
        stopPreview();
        if (gd.wasCanceled()) {
            if (imp != null) {
                imp.setRoi(originalRoi);
            }
            valueTolerance = saveValueTolerance;
            gradientTolerance = saveGradientTolerance;
            colorSensitivity = saveColorSensitivity;
            connection = saveConnection;
            includeHoles = saveIncludeHoles;
            useEyedropper = saveUseEyedropper;
        }
    }

    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        valueTolerance = gd.getNextNumber();
        colorSensitivity = 0.01 * gd.getNextNumber();
        gradientTolerance = gd.getNextNumber();
        useEyedropper = gd.getNextBoolean();
        connection = gd.getNextChoiceIndex();
        includeHoles = gd.getNextBoolean();
        if (gd.invalidNumber() || valueTolerance < 0 || colorSensitivity < -1 || colorSensitivity > 1 || gradientTolerance < 0) {
            return false;
        }
        if (previewing) {
            synchronized (this) {
                previewThread.interrupt(); //trigger preview
            }
        }
        return true;
    }

    /**
     * The tool icon; for details see
     * http://rsb.info.nih.gov/ij/developer/macro/macros.html#icons
     */
    public String getToolIcon() {
        return "C333F8082C555F8282C777F8482C999F8682CbbbF8882CdddF8a82Cf00Lee44O1133";
    }

    /**
     * This method can be called from a macro or plugin to do the wand operation
     * in the current thread. In contrast to the operation triggered by the
     * mode, this does not lock the image. 'options' can include a Strong for
     * the connection mode ("4-connected", "non-contiguous"; "8-connected" is
     * default and 'eydropper' for using the foreground color (eydropper color),
     * and "include" for including holes).
     */
    public static void doWand(int xStart, int yStart, int tolerence, String options) {
        ImagePlus imp = WindowManager.getCurrentImage();
        (new Contour_Roi_Tool()).doWand(imp, xStart, yStart, tolerence, options);
    }

    private void doWand(ImagePlus imp, int xStart, int yStart, double tolerence, String options) {
        if (imp == null) {
            return;
        }
        valueTolerance = tolerence;
        if (options != null) {
            options = options.toLowerCase();
            for (int i = 0; i < CONNECTED_STRINGS.length; i++) {
                if (options.contains(CONNECTED_STRINGS[i])) {
                    connection = i;     //found a keyword for connection mode
                }
            }
            useEyedropper = options.contains("eydropper");
            includeHoles = options.contains("include");
        }
        imp.saveRoi();
        doWand(imp, xStart, yStart);
    }

    /**
     * The wand operation, triggered by the mouse. Since the wand operation may
     * take a while, it is performed in a separate thread (not the event queue),
     * so that the progress bar still works and ImageJ does not become
     * unresponsive. The image is locked during the operation.
     */
    public void mousePressed(ImagePlus imp, MouseEvent e) {
        ImageCanvas ic = imp.getCanvas();
        xStart = ic.offScreenX(e.getX());
        yStart = ic.offScreenY(e.getY());
        yshift = 0;
        valueTolerance0 = 255 * (imp.getProcessor().getPixelValue(xStart, yStart) / imp.getStatistics().max);
        valueTolerance = valueTolerance0;
        update(imp);
    }

    public void mouseDragged(ImagePlus imp, MouseEvent e) {
        ImageCanvas ic = imp.getCanvas();
        yshift = yStart - ic.offScreenY(e.getY());
        valueTolerance = valueTolerance0 + yshift;
        valueTolerance = Math.min(254, Math.max(1, valueTolerance));
        update(imp);
    }

    private void update(final ImagePlus imp) {
        if (imp.isLocked()) {
            return; //image has been locked previously
        } else {
            imp.lock();
        }
        imp.saveRoi();
        final int x = xStart, y = yStart;
        Thread doWandThread = new Thread(new Runnable() {
            public void run() {
                doWand(imp, x, y);
                imp.unlock();
            }
        });
        doWandThread.setName("Versatile Wand");
        doWandThread.start();

        if (Recorder.record) {
            String options = CONNECTED_STRINGS[connection];
            if (useEyedropper) {
                options += " eydropper";
            }
            if (includeHoles) {
                options += " include";
            }
            if (Recorder.scriptMode()) {
                Recorder.recordCall("vwt = IJ.runPlugIn(\"Contour_Roi_Tool\", \"\");");
                Recorder.recordCall("vwt.doWand(imp, " + xStart + ", " + yStart + ", " + valueTolerance + ", \"" + options + "\");");
            } else {
                String s = "call(\"Contour_Roi_Tool.doWand\", " + xStart + ", " + yStart + ", " + (int) valueTolerance + ", \"" + options + "\");\n";
                Recorder.recordString(s);
            }
        }
    }

    //int iPxlForDebug;
    /**
     * Here the wand operation actually happens
     */
    private void doWand(ImagePlus imp, int x0, int y0) {
        //try{
        //IJ.log("doWand at "+x0+","+y0);
        IJ.showProgress(0.0);
        boolean shiftKeyDown = IJ.shiftKeyDown();
        boolean altKeyDown = IJ.altKeyDown();
        ImageProcessor ip = imp.getProcessor();
        int width = ip.getWidth();
        int height = ip.getHeight();

        // prepare mask
        ByteProcessor maskIp = new ByteProcessor(width, height);
        byte[] mPixels = (byte[]) maskIp.getPixels();
        // offsets to neighboring pixels in x, y and pixel number
        int[] dirXoffset = new int[]{0, 1, 1, 1, 0, -1, -1, -1};
        int[] dirYoffset = new int[]{-1, -1, 0, 1, 1, 1, 0, -1,};
        int[] dirOffset = new int[]{-width, -width + 1, +1, +width + 1, +width, +width - 1, -1, -width - 1};
        // pixelPointers array for positions that we have to process; initial size is 4096 = 0x1000 points
        int pixelPointerMask = 0xf;   //binary AND with this mask gets array index
        int[] pixelPointers = new int[pixelPointerMask + 1];

        // Pixel value range
        Calibration cal = imp.getCalibration();
        float lowLimit, highLimit;
        float grayRef;
        Color eyeDropperColor = null;
        if (useEyedropper) {
            eyeDropperColor = Toolbar.getForegroundColor(); // reference color is eyedropper color
            int index = ip.getBestIndex(eyeDropperColor);
            float[] calTable = ip.getCalibrationTable();
            grayRef = calTable == null ? index : calTable[index];
        } else {
            grayRef = ip.getPixelValue(x0, y0);             // reference color from this pixel
        }
        if (ip.getMinThreshold() == ImageProcessor.NO_THRESHOLD) {
            lowLimit = grayRef - (float) valueTolerance;
            //highLimit = grayRef + (float) valueTolerance;
            highLimit = (float) ip.getStatistics().max;
        } else if (cal == null) {
            lowLimit = (float) ip.getMinThreshold();
            highLimit = (float) ip.getMaxThreshold();
        } else {
            lowLimit = (float) cal.getCValue((int) (ip.getMinThreshold() + 0.5));
            highLimit = (float) cal.getCValue((int) (ip.getMaxThreshold() + 0.5));
        }
        if (grayRef < lowLimit || grayRef > highLimit) {
            IJ.showStatus("Not in Thresholded Area");
            IJ.beep();
            return;
        }

        // Prepare color threshold
        boolean colorMode = (ip instanceof ColorProcessor) && colorSensitivity > -1;
        int[] rgbPixels = null;
        int r0 = 0, g0 = 0, b0 = 0;                       //rgb of the starting/reference point
        double rWeight = 0, gWeight = 0, bWeight = 0;     //weights for getting 'parallel' or 'gray' comonent from rgb
        if (colorMode) {
            rgbPixels = (int[]) ip.getPixels();
            if (useEyedropper) {
                float[] rgbalpha = eyeDropperColor.getRGBComponents(null);
                r0 = (int) (rgbalpha[0] * 255 + 0.5);
                g0 = (int) (rgbalpha[1] * 255 + 0.5);
                b0 = (int) (rgbalpha[2] * 255 + 0.5);
            } else {
                int c0 = rgbPixels[x0 + y0 * width];
                r0 = (c0 >> 16) & 0xff;
                g0 = (c0 >> 8) & 0xff;
                b0 = (c0) & 0xff;
            }
            double[] refColor = ((r0 == 0 && g0 == 0 && b0 == 0) || colorSensitivity <= 0)
                    ? ((ColorProcessor) ip).getWeightingFactors() : new double[]{r0, g0, b0};
            normalize(refColor);
            rWeight = refColor[0];
            gWeight = refColor[1];
            bWeight = refColor[2];
        }

        // Prepare gradient threshold
        boolean useGradient = gradientTolerance > 0 && gradientTolerance < valueTolerance && connection != NON_CONTIGUOUS;
        float aspectRatioSqr = 1;   // (pixel width/pixel height)^2
        float toleranceGrayGradTmp = (float) gradientTolerance;
        if (cal != null) {
            toleranceGrayGradTmp *= (float) cal.pixelWidth;
            aspectRatioSqr = (float) (cal.pixelWidth / cal.pixelHeight);
            aspectRatioSqr *= aspectRatioSqr;
        }
        float toleranceGrayGrad2 = toleranceGrayGradTmp * toleranceGrayGradTmp;

        int ymin = height;          // the uppermost selected pixel, needed for reliable IJ.doWand
        if (connection == NON_CONTIGUOUS) {
            //create mask from full image

            if (colorMode) {
                for (int i = 0; i < width * height; i++) //{iPxlForDebug=i;
                {
                    if (checkColor(rgbPixels[i], r0, g0, b0, rWeight, gWeight, bWeight)) {
                        mPixels[i] = INSIDE;
                    }
                }
                //}
            } else {
                for (int y = 0, offset = 0; y < height; y++) {
                    for (int x = 0; x < width; x++, offset++) {
                        float v = ip.getPixelValue(x, y);
                        if (v >= lowLimit && v <= highLimit) {
                            mPixels[offset] = INSIDE;
                        }
                    }
                }
            }
        } else {
            //queue-based flood fill algorithm
            int lastCoord = 0;
            int offset0 = x0 + y0 * width;
            mPixels[offset0] = INSIDE;
            pixelPointers[0] = offset0;
            for (int iCoord = 0; iCoord <= lastCoord; iCoord++) {
                int offset = pixelPointers[iCoord & pixelPointerMask];
                int x = offset % width;
                int y = offset / width;
                boolean isInner = (x != 0 && y != 0 && x != (width - 1) && y != (height - 1));
                float v = ip.getPixelValue(x, y);
                boolean largeGradient = false;
                float xGradient = 0, yGradient = 0;
                if (useGradient) {
                    if (isInner) {
                        float vpp = ip.getPixelValue(x + 1, y + 1);
                        float vpm = ip.getPixelValue(x + 1, y - 1);
                        float vmp = ip.getPixelValue(x - 1, y + 1);
                        float vmm = ip.getPixelValue(x - 1, y - 1);

                        xGradient = 0.125f * ( //Sobel-filter like gradient
                                2f * (ip.getPixelValue(x + 1, y) - ip.getPixelValue(x - 1, y))
                                + vpp - vmm + (vpm - vmp)); // v(x+1) - v(x-1)
                        yGradient = 0.125f * (2f * (ip.getPixelValue(x, y + 1) - ip.getPixelValue(x, y - 1))
                                + vpp - vmm - (vpm - vmp)); // v(y+1) - v(y-1)
                    } else {
                        int xCount = 0, yCount = 0;
                        for (int d = 0; d < 8; d++) {
                            if (isWithin(ip, x, y, d)) {
                                int x2 = x + dirXoffset[d];
                                int y2 = y + dirYoffset[d];
                                float v2 = ip.getPixelValue(x2, y2);
                                int weight = (2 - (d & 0x1));   //2 for straight, 1 for diag
                                xGradient += dirXoffset[d] * (v2 - v) * weight;
                                xCount += weight * (dirXoffset[d] != 0 ? 1 : 0);
                                yGradient += dirYoffset[d] * (v2 - v) * weight;
                                yCount += weight * (dirYoffset[d] != 0 ? 1 : 0);
                            }
                        }
                        xGradient /= xCount;
                        yGradient /= yCount;
                    }
                    largeGradient = xGradient * xGradient + yGradient * yGradient * aspectRatioSqr > toleranceGrayGrad2;
                }
                for (int d = 0; d < 8; d += (connection == FOUR_CONNECTED) ? 2 : 1) {  //analyze all neighbors (in 4 or 8 directions)
                    int offset2 = offset + dirOffset[d];
                    if ((isInner || isWithin(ip, x, y, d)) && mPixels[offset2] == UNKNOWN) {
                        int x2 = x + dirXoffset[d];
                        int y2 = y + dirYoffset[d];
                        float v2 = 0;
                        boolean valueOK;
                        if (largeGradient || !colorMode) {
                            v2 = ip.getPixelValue(x2, y2);
                        }
                        if (colorMode) {
                            int c2 = rgbPixels[offset2];
                            valueOK = checkColor(c2, r0, g0, b0, rWeight, gWeight, bWeight);
                        } else // if !colorMode
                        {
                            valueOK = v2 >= lowLimit && v2 <= highLimit;
                        }
                        if (!valueOK) {
                            mPixels[offset2] = OUTSIDE;                 //don't analyze any more
                        } else if (!largeGradient || (v2 - v) * (xGradient * dirXoffset[d] + yGradient * dirYoffset[d]) <= 0) {
                            mPixels[offset2] = INSIDE;                  //add new point
                            if (ymin > y2) {
                                ymin = y2;
                            }
                            if (lastCoord - iCoord > pixelPointerMask) {  //not enough space in array, expand it
                                int newSize = 2 * (pixelPointerMask + 1);
                                int newMask = newSize - 1;
                                int[] newPixelPointers = new int[newSize];
                                System.arraycopy(pixelPointers, 0, newPixelPointers, 0, pixelPointerMask + 1);
                                System.arraycopy(pixelPointers, 0, newPixelPointers, pixelPointerMask + 1, pixelPointerMask + 1);
                                pixelPointers = newPixelPointers;
                                pixelPointerMask = newMask;
                            }
                            lastCoord++;
                            pixelPointers[lastCoord & pixelPointerMask] = offset2;

                        }
                    }
                } //for direction d
                if ((iCoord & 0xfff) == 1) {
                    IJ.showProgress(0.5 * Math.sqrt(iCoord / (double) (width * height))); //faster progress bar at the beginning
                    if (Thread.currentThread().isInterrupted()) {
                        return;
                    }
                }
            } //for iCoord
        } //else: if contiguous

        if (Thread.currentThread().isInterrupted()) {
            return;
        }
        IJ.showProgress(0.5);
        Roi previousRoi = imp.getRoi();
        //convert mask to roi
        Roi roi = null;
        if (includeHoles && connection != NON_CONTIGUOUS) {
            Wand w = new Wand(maskIp);
            w.autoOutline(0, ymin, 255, 255);   //starting at x0, y0 might trace an inner hole
            if (w.npoints > 0) {
                roi = new PolygonRoi(w.xpoints, w.ypoints, w.npoints, Roi.TRACED_ROI);
                if (Thread.currentThread().isInterrupted()) {
                    return;
                }
                imp.killRoi();
                imp.setRoi(roi);
            }
        } else {
            maskIp.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE);
            if (Thread.currentThread().isInterrupted()) {
                return;
            }
            imp.killRoi();
            ThresholdToSelection tts = new ThresholdToSelection();
            tts.setup("", imp);
            tts.run(maskIp);
            IJ.showStatus("");
            roi = imp.getRoi();
        }

        // add/subtract this ROI to the previous one if the shift/alt key is down
        if (Thread.currentThread().isInterrupted()) {
            return;
        }
        if (roi != null && previousRoi != null) {
            roi.update(shiftKeyDown, altKeyDown);
        }

        IJ.showProgress(1.0);
        //}catch(Exception e){IJ.handleException(e);}
    }

    /**
     * Returns whether pixel with rgb is inside color tolerance
     */
    boolean checkColor(int rgb, int r0, int g0, int b0, double rWeight, double gWeight, double bWeight) {
        int deltaR = ((rgb >> 16) & 0xff) - r0;
        int deltaG = ((rgb >> 8) & 0xff) - g0;
        int deltaB = ((rgb) & 0xff) - b0;
        double deltaSqr = (deltaR * deltaR + deltaG * deltaG + deltaB * deltaB);
        //if(iPxlForDebug==253740)IJ.log("cSens="+(float)colorSensitivity+" dR="+deltaR+" dG="+deltaG+" dB="+deltaB);
        //if(iPxlForDebug==253740)IJ.log("rWeight="+(float)rWeight+" gWeight="+(float)gWeight+" bWeight="+(float)bWeight);
        if (colorSensitivity == 0) {
            return (1. / 3.) * deltaSqr <= valueTolerance * valueTolerance;
        } else {
            double deltaPar = deltaR * rWeight + deltaG * gWeight + deltaB * bWeight;
            double deltaParSqr = deltaPar * deltaPar;
            // we want the following weighting of delta-squared parallel and perpendicular components:
            //   parallel + (1+colorSensitivity)*perpendicular if colorSensitivity <0,
            //   parallel*(1-colorSensitivity) + perpendicular if colorSensitivity >0
            // The perpendicular component is deltaSqr - deltaParSqr.
            //if(iPxlForDebug==253740)IJ.log("deltaSqr="+(float)deltaSqr+" deltaParSqr="+(float)deltaParSqr+" deltaPerpSqr="+(float)(deltaSqr-deltaParSqr));
            if (colorSensitivity < 0) {
                return deltaSqr * (1. + colorSensitivity) - deltaParSqr * colorSensitivity <= valueTolerance * valueTolerance;
            } else {
                return deltaSqr - deltaParSqr * colorSensitivity <= valueTolerance * valueTolerance;
            }
        }
    }

    /**
     * Returns whether the neighbor in a given direction is within the image
     * NOTE: it is assumed that the pixel x,y itself is within the image! Uses
     * class variables width, height: dimensions of the image
     *
     * @param x x-coordinate of the pixel that has a neighbor in the given
     * direction
     * @param y y-coordinate of the pixel that has a neighbor in the given
     * direction
     * @param direction the direction from the pixel towards the neighbor (see
     * makeDirectionOffsets)
     * @return true if the neighbor is within the image (provided that x, y is
     * within)
     */
    boolean isWithin(ImageProcessor ip, int x, int y, int direction) {
        int width = ip.getWidth();
        int height = ip.getHeight();
        int xmax = width - 1;
        int ymax = height - 1;
        switch (direction) {
            case 0:
                return (y > 0);
            case 1:
                return (x < xmax && y > 0);
            case 2:
                return (x < xmax);
            case 3:
                return (x < xmax && y < ymax);
            case 4:
                return (y < ymax);
            case 5:
                return (x > 0 && y < ymax);
            case 6:
                return (x > 0);
            case 7:
                return (x > 0 && y > 0);
        }
        return false;   //should never occur, we use it only for directions 0-7
    } // isWithin

    /* normalize 3-element vector */
    void normalize(double[] v) {
        double norm = 1. / Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        for (int i = 0; i < 3; i++) {
            v[i] *= norm;
        }
        //IJ.log("weights:"+(float)v[0]+","+(float)v[1]+","+(float)v[2]);
    }

    /**
     * Runs the preview thread while the dialog is open.
     */
    private void startPreview(final ImagePlus imp, final int xStart, final int yStart) {
        previewThread = new Thread(new Runnable() {
            public void run() {
                previewing = true;
                do {
                    doWand(imp, xStart, yStart);
                    synchronized (this) {
                        if (!previewing) {
                            return;
                        }
                        if (!Thread.interrupted()) {
                            try {
                                wait();
                            } catch (InterruptedException e) {
                            }
                        }
                    }
                } while (previewing);
            }
        });
        previewThread.setName("Versatile Wand Preview");
        previewThread.start();
    }

    private void stopPreview() {
        synchronized (this) {
            previewing = false;
            if (previewThread != null) {
                previewThread.interrupt();
            }
            IJ.showProgress(1.0);
        }
    }

}
