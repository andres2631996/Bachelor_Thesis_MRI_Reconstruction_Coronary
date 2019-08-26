Function RegionGrowing, Coord, Img, threshold
; Segmentation method based on introducing an initial coordinate (coord), supposed to be in the object
; It determines whether the surrounding pixels belong to the object, too, seeing if their values are inside a threshold
; In each iteration, more pixels are evaluated and if they belong to the object, it "grows"
; As a consequence, all the pixels in the foreground are spatially connected between them
; Inputs:
; Coord: initial coordinate where the segmentation starts
; Img: image to be segmented by region growing
; Threshold: vector with the minimum and maximum values that are allowed to make part of the object
; Output:
; Mask: segmentation mask of the size of the input image with the object
outputRes=size(Img, /dimensions) ; Matrix size
mask=bytarr(outputRes(0),outputRes(1)) ; Final segmentation mask
; EXCEPTIONS
if (Img(Coord(0),Coord(1)) lt threshold(0))||(Img(Coord(0),Coord(1)) gt threshold(1)) then begin
	mask=make_array(outputRes(0),outputRes(1),value=0)
	print, mask ; This happens if we select an initial coordinate lying outside the object to be segmented
endif else begin
    mask(Coord(0),Coord(1))=1
endelse
if (Coord(0) lt 0) || (Coord(1) lt 0) || (Coord(0) gt outputRes(0)-1) || (Coord(1) gt outputRes(1)-1) then begin
 	mask=make_array(outputRes(0),outputRes(1),value=0)
	print, mask ; This happens if we select a coordinate which is negative or greater than the matrix size
endif
; SEGMENTATION METHOD
cont=fix(1) ; Iteration counter
res=bytarr((2*cont+1)^2) ; Array with the segmentation results in each iteration
cont_element=fix(0) ; Element counter
while (total(res) gt 0) || (total(res) lt outputRes(0)*outputRes(1)) do begin ; The loop ends when we do not find any more pixels to add to the object or when all the image belongs to the object
; (This last condition is stated in order not to have an infinite loop, just in case)
	for i=-cont,cont do begin ; Explore rows
		for j=-cont,cont do begin ; Explore columns
			if (Coord(0)+i gt outputRes(0)-1) || (Coord(1)+j gt outputRes(1)-1) || (Coord(0)+i lt 0) || (Coord(1)+j lt 0) then begin
				res(cont_element)=0 ; This happens if the object is touching the edges of the image, so the supposed coordinates of the mask exceed the matrix size or are negative
			endif else begin ; In this case, the coordinates where we are looking for are inside the image coordinates
		 		if (Img(Coord(0)+i, Coord(1)+j) ge threshold(0)) && (Img(Coord(0)+i, Coord(1)+j) le threshold(1)) then begin
		 			mask(Coord(0)+i, Coord(1)+j)=1
		 			res(cont_element)=1
		 		endif else begin
	                mask(Coord(0)+i, Coord(1)+j)=0
	                res(cont_element)=0
	            endelse
	        endelse
	        cont_element=cont_element+1
	    endfor
	endfor
	cont_element=0
	cont=cont+1
	res=bytarr((2*cont+1)^2)
endwhile
return, mask
end