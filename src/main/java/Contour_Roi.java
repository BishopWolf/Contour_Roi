
import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.tool.PlugInTool;
import ij.plugin.filter.ThresholdToSelection;
import java.awt.*;
import java.awt.event.*;

/**
 * An ImageJ contour Roi tool.
 *
 * @author alex.vergara
 */
public class Contour_Roi extends PlugInTool {

    private double max = 0;
    private Roi finalRoi = null;
    private boolean MouseDownState = false;
    // remember previous position for level adjust
    private int xStart = -1, yStart = -1;

    /**
     * The tool icon; for details see
     * http://rsb.info.nih.gov/ij/developer/macro/macros.html#icons
     */
    @Override
    public String getToolIcon() {
        return "C333F8082C555F8282C777F8482C999F8682CbbbF8882CdddF8a82Cf00Lee44O1133";
    }

    /**
     * This method can be called from a macro or plugin to do the contour
     * operation in the current thread.
     *
     * @param xStart current mouse x position
     * @param yStart current mouse y position
     * @param level current level
     */
    public static void doContour(int xStart, int yStart, int level) {
        ImagePlus imp = WindowManager.getCurrentImage();
        (new Contour_Roi()).doContour(imp, xStart, yStart, level);
    }

    private void doContour(ImagePlus imp, int xStart, int yStart, int level) {
        if (imp == null) {
            return;
        }
        Point p = new Point(xStart, yStart);
        this.finalRoi = getObject(imp, p, level, this.max);
        if (this.finalRoi != null) {
            imp.deleteRoi();
            imp.setRoi(this.finalRoi);
        }
    }

    @Override
    public void mousePressed(ImagePlus imp, MouseEvent e) {
        this.MouseDownState = true;
        this.max = imp.getStatistics().max;
    }

    @Override
    public void mouseReleased(ImagePlus imp, MouseEvent e) {
        this.MouseDownState = false;
    }

    /**
     * The Roi operation, triggered by the mouse move.
     *
     * @param imp
     * @param e
     */
    @Override
    public void mouseMoved(final ImagePlus imp, MouseEvent e) {
        if (!imp.lock()) {
            return; //image has been locked previously
        }
        if (!MouseDownState) {
            return; // mouse not in down state
        }
        Thread.currentThread().interrupt();
        //imp.saveRoi();
        ImageCanvas ic = imp.getCanvas();
        xStart = ic.offScreenX(e.getX());
        yStart = ic.offScreenY(e.getY());
        final int x = xStart, y = yStart;
        final int level = (int) imp.getProcessor().getPixelValue(x, y);
        Thread doRoiThread = new Thread(new Runnable() {
            public void run() {
                doContour(imp, x, y, level);
                imp.unlock();
            }
        });
        doRoiThread.setName("Contour Roi");
        doRoiThread.start();

    }

    /**
     *
     * @param imp The image object
     * @param p the point on which we calculate the boundary, this point is
     * always inside final roi
     * @param level the level to get isocontour
     * @param max the maximum of the image
     * @return the isocontour roi at specified level containing the input point
     */
    public static Roi getObject(ImagePlus imp, Point p, double level, double max) {
        ImageProcessor ip2 = imp.getProcessor().duplicate();
        ip2.setThreshold(level, max, ImageProcessor.BLACK_AND_WHITE_LUT);
        ThresholdToSelection ts = new ThresholdToSelection();
        Roi roi = ts.convert(ip2);
        ShapeRoi Sroi = new ShapeRoi(roi);
        Roi[] listroi = Sroi.getRois();
        for (Roi r : listroi) {
            if (Thread.currentThread().isInterrupted()) {
                return null;
            }
            if (r.contains((int) p.getX(), (int) p.getY())) {
                roi = r;
                break;
            }
        }
        return roi;
    }

}
