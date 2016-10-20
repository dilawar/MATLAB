function multiROI_CROIEditor
myimage = imread('/Users/ananthamurthy/Desktop/Work/Imaging/20131218/cells-spnt2-ROI-1r.tif',1);
myimage=imadjust(myimage);
roiwindow = CROIEditor(myimage);

addlistener(roiwindow,'MaskDefined',@your_roi_defined_callback)

    function your_roi_defined_callback(h,e)
        [mask, labels, n] = roiwindow.getROIData;
        delete(roiwindow);
    end
end