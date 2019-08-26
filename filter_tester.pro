function filter_tester, outputRes, mask
maskKspace=complexarr(outputRes(0), outputRes(1))
maskKspace=mask(*,*,0)
maskKspace=Generate_KSpace_mask(maskKspace) ; Obtain a mask to select where in the k-space we are going to work
vectorX=(where(maskKspace(*,0) eq 0)) ; Indexes where the rows of the mask are 0
vectorY=(where(maskKspace(0,*) eq 0)) ; Indexes where the columns of the mask are 0
radio=float([vectorX(0),vectorY(0)])+1 ; Contains the first indexes of rows and columns with zeros. A one is added
LowerRadio=radio-[3.,3.] ; The lower radio is 3 elements shorter than the initial radio
HalfScan=[vectorX(n_elements(vectorX)-1),vectorY(n_elements(vectorY)-1)] ; Halfscan contains the last elements of the rows and columns of the mask which are 0
LowerHalfScan=HalfScan+[3.,3.] ; Lower half scan is 3 elements higher than the original half scan
maskKspace=fltarr(outputRes(0),outputRes(1)) ; Float array with dimensions the matrix size of the final image
filter=fltarr(outputRes(0),outputRes(1)) ; Float array with dimensions the matrix size of the final image
for j=0, outputRes(1)-1 do begin
   for i=0, outputRes(0)-1 do begin
        radius=((i-outputRes(0)/2.)/radio(0))^2.+((j-outputRes(1)/2.)/radio(1))^2. ; Distance to the central pixel in the k-space normalized by the first elements where the mask is 0
        if radius le 1 then maskKspace(i,j)=1.0 ; Looks if the radius is lower than 1, assigning then the mask a 1 --> LP filter in the k-space
        LowerRadius=((i-outputRes(0)/2.)/LowerRadio(0))^2.+((j-outputRes(1)/2.)/LowerRadio(1))^2. ; Distance to the central pixel in the k-space normalized by the first elements where the mask is 0, minus 3
        if LowerRadius le 1 then filter(i,j)=1.0 ; Looks if the lower radius is lower than 1, assigning then the filter a 1 --> LP filter in the k-space
    endfor
endfor

maskKspace=shift(temporary(maskKspace),outputRes(0)/2.,outputRes(1)/2.) ; Shift the mask with half of the dimensions of the matrix size
maskKspace(radio(0):HalfScan(0),*)=0 ; The mask between the first and last indexes of rows and columns where it was 0 is now 0
maskKspace(*,radio(1):HalfScan(1))=0

filter=shift(temporary(filter),outputRes(0)/2.,outputRes(1)/2.) ; Shift the filter with half of the dimensions of the matrix size
filter(LowerRadio(0):LowerHalfScan(0),*)=0 ; The filter between the first indexes of rows and columns where it was 0, minus 3, and the last indexes of rows and columns where it was 0, plus 3, is now 0
filter(*,LowerRadio(1):LowerHalfScan(1))=0

HomodineFilter=maskKspace ; The final mask is a LP Homodyne filter

filter=shift(temporary(filter),-outputRes(0)/2.,-outputRes(1)/2.) ; The filter is unshifted
filter=(convol(temporary(filter), FilterGen(25,19),/EDGE_t)) ; The filter is now convolved with a Gaussian filter of 25 samples and a width of 19, convolving also in the edges

maskKspace=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.) ; The mask is unshifted

weight=fltarr(1,outputRes(1)) ; Vector of weights with the dimensions of Z
weight(0,*)=exp(-(findgen(outputRes(1))-outputRes(1)/2.)^2./40000) ; Weights are normally distributed
weight=rebin(temporary(weight), outputRes(0), outputRes(1)) ; Weights are now distributed in a matrix with size the matrix size of our image

a=fltarr(outputRes(0),1) ; Vector of weights with the dimensions of Y
a(*,0)=exp(-(findgen(outputRes(0))-outputRes(0)/2.)^2./40000) ; Weights are normally distributed
a=rebin(a, outputRes(0), outputRes(1)) ; Weights are now distributed in a matrix with size the matrix size of our image

weight = temporary(weight)*a ; Weights in Y and Z directions are now multiplied between them

filter=temporary(filter)*shift(HomodineFilter, -outputRes(0)/2.,-outputRes(1)/2.)*weight
; Homodyne filter (the initial filter) is shifted with half of the dimensions of the final matrix size and multiplied by the convolved version of the filter and by the Gaussian weights
weight=0

filter=reform(shift(temporary(filter), -outputRes(0)/2.,-outputRes(1)/2.),  outputRes(0),outputRes(1),1,1) ; The filter is shifted again and adapted to the matrix size of the final image
filter=rebin(temporary(filter),outputRes(0),outputRes(1),1,outputRes(3)) ; The filter is adapted to the dimensions of the final matrix size and the slice number

return, filter
end