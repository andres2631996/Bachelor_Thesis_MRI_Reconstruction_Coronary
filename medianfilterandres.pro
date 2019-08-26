Function MedianFilterAndres, outputRes, volumen, XYwidth
; Function which computes the median filtering of an image by introducing the image to be filtered and the size of the kernel
; Image is a 2D image with a certain matrix size (M x N) and kernel is a square matrix of odd dimensions
NewImg=fltarr(outputRes(0)+XYwidth(0)-1, outputRes(1)+XYwidth(1)-1) ; Auxiliary image used to filter the edges
volumenFilter=fltarr(outputRes(0), outputRes(1), outputRes(3))
for k=0, outputRes(3)-1 do begin ; Go through all slices
	NewImg(floor(XYwidth(0)/2):outputRes(0)+floor(XYwidth(0)/2)-1, floor(XYwidth(1)/2):outputRes(1)+floor(XYwidth(1)/2)-1)=abs(volumen(*,*,k))
	for i=0, outputRes(1)-1 do begin ; Go through all columns
		NewImg(0:floor(XYwidth(0)/2)-1,i+floor(XYwidth(1)/2)-1)=make_array(floor(XYwidth(0)/2),value=abs(volumen(0,i,k))) ; Values to the left of the edge in the original image
		NewImg(outputRes(0):outputRes(0)+floor(XYwidth(0)/2)-1,i+floor(XYwidth(1)/2)-1)=make_array(floor(XYwidth(0)/2), value=abs(volumen(outputRes(0)-1,i,k))) ; Values to the right of the edge in the original image
	endfor
	for j=0, outputRes(0)-1 do begin ; Go through all rows
		NewImg(j+floor(XYwidth(0)/2)-1, 0:floor(XYwidth(1)/2)-1)=make_array(1, floor(XYwidth(1)/2),value=abs(volumen(j,0,k))) ; Values over the edge in the original image
		NewImg(j+floor(XYwidth(0)/2)-1, outputRes(1):outputRes(1)+floor(XYwidth(1)/2)-1)=make_array(1, floor(XYwidth(1)/2), value=abs(volumen(j,outputRes(1)-1,k))) ; Values below the edge in the original image
	endfor
	for i=0, outputRes(0)-1 do begin
		for j=0, outputRes(1)-1 do begin
			volumenFilter(i,j,k)=median(reform(NewImg(i:i+XYwidth(0)-1,j:j+XYwidth(1)-1),XYwidth(0)*XYwidth(1),1))
		endfor
	endfor
endfor
return, volumenFilter
end