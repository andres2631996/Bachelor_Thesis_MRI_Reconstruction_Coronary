pro whole_tester
; Receives as input the header of the cpx, the data, outputRes, the slices, the information of the Scan protocol, the SENSE factor, the mask and the time mask
; DATA OBTENTION
dir='H:\CORONARIAS\' ; Directory where to read the raw data.
dirfinal=dir+'FinalVolume\'
DataHeader=file_search(dir+'raw*.list')
DataData=file_search(dir+'raw*.data')
SlicesCenter = 0
CoilHeaderCPX=file_search(dir+'cpx*.list')
CoilDataCPX=file_search(dir+'cpx*.data')
NCoils=16
Coils=read_cpx(CoilHeaderCPX, CoilDataCPX, QBodyCoil, NCoils)
outputRes=float([324,196, NCoils, 14])
RefProto=ReadProtocol(dir+'Ref_XL_Torso.txt')
ScanProto=ReadProtocol(dir+'CORONARIAS+EDEMA.txt') ; Same for the Scan protocol
FOVShift=(ceil(ScanProto.SENSE)-ScanProto.SENSE)/ceil(ScanProto.SENSE) ; Shift of the FOV: depends on the SENSE factor
FOVCorrection=scanproto.FOV ; Final FOV
FOVCorrection=scanproto.FOV/FOVCorrection ; Correction of FOV
ScanProto.SENSE=ceil(ScanProto.SENSE) ; The SENSE factor is approximated to the closer and higher integer
SENSEFactor=ScanProto.SENSE
QBThreshold=0.05
; DATA AND MASK OBTENTION
mask=mask_tester(DataHeader, DataData, outputRes,ScanProto, SENSEFactor)
data=data_tester(DataHeader, DataData, outputRes, SlicesCenter, ScanProto, ScanProto.SENSE, mask)

DataHeader=0
DataData=0
CoilHeaderCPX=0
CoilDataCPX=0
; SENSITIVITY MAPS, Q-BODY COIL AND EXPORT MASK OBTENTION
QBodyCoil=QBodyCoilTester(Coils, RefProto, ScanProto, outputRes, FOVCorrection, QBodyCoil, SlicesCenter)
ExportMask=ExportMask_tester(QBodyCoil, QBThreshold, outputRes)
SenseMaps=SenseMapsTester(Coils, QBodyCoil, RefProto, ScanProto, outputRes, FOVCorrection, SlicesCenter, mask, ExportMask)
QBodyCoil=temporary(QBodyCoil)*ExportMask ; The Q-body coil is also masked by the export mask
Coils=0
; FILTERING
filter=filter_tester(outputRes, mask)
for i=0, outputRes(2)-1 do begin
    data(*,*,i,*)=temporary(data(*,*,i,*))*filter ; Filtered data
endfor
filter=reform(filter) ; Dimensions with just one value are removed
mask=fltarr(outputRes(0),outputRes(1),outputRes(2)) ; Float array analyzed in just one slice
ind=where(abs(data(*,*,*,0)), count) ; Indexes of the absolute value of our data
mask(ind)=1. ; The mask is one in all the indexes of our data
;mask=total(mask,5) ; Sum of all the elements of data for the same cardiac phase
ind=where(abs(mask), count) ; Keep the mask indexes

; DATA SAMPLING
Data=temporary(Data(ScanProto.SENSE(0)*findgen(outputRes(0)/(ScanProto.SENSE(0))),*,*,*)) ; Data is restricted in Y and Z directions by the SENSE factors in both directions
Data=temporary(Data(*,ScanProto.SENSE(1)*findgen(outputRes(1)/(ScanProto.SENSE(1))),*,*))

; COIL SELECTION
;ChannelSNR=total(total(abs(data(0,0,*,*)),2),2)/outputRes(3) ; Sum of the absolute value of the first element (Row 0, Column 0) in all the slices, divided by the number of slices and multiplied by the number of cardiac phases
;ind=reverse(sort(ChannelSNR)) ; Absolute values of the first elements of data in each slice sorted in descending order
;ChannelSelection=where(ChannelSNR gt (mean(ChannelSNR(ind))-2.0*stddev(ChannelSNR(ind))), NumSelectedChannels)
;; Indexes where the absolute values of the first elements of the data in each slice are higher than the sorted absolute values of the first elements (from the first slice to the slice 2/3 of the of the total number of coils)
;; minus twice the standard deviation of the sorted absolute values (from the first element to the slice whose number is equal to half of the coils used
;; The number of times this happens is saved in "NumSelectedChannels"
;if NumSelectedChannels lt 8 then NumSelectedChannels=outputRes(2)
;; See if NumSelectedChannels is lower than 8. In that case, the number of selected channels equals the number of coils
;if NumSelectedChannels lt outputRes(2) then begin
;	;print, ChannelSelection ; Print the channel indexes if the number of selected channels is lower than the number of coils
;	ChannelSNR=0
;	data=data(*,*,ChannelSelection,*) ; The data is now just kept in the channels we selected
;	SenseMaps=temporary(SenseMaps(*,*,ChannelSelection,*)) ; Sensitivity Maps are saved temporarily in the selected channels
;    outputRes(2)=NumSelectedChannels ; The number of channels in outputRes is equal to the number of selected channels
;endif
ScanProto=0
RefProto=0
; FIRST RECONSTRUCTION
SVDthreshold=1e-1
volumen1=final_reconstruction(SVDthreshold, outputRes, data, SENSEFactor, SenseMaps, QBodyCoil)
volumen1=temporary(volumen1)*abs(ExportMask)
; MEDIAN FILTERING
XYwidth=intarr(3)
XYwidth=[5,5,5]
QBodyCoil1=complexarr(outputRes(0),outputRes(1),outputRes(3))
QBodyCoil1=QBodyCoil
QBodyCoil=MedianFilterAndres(abs(volumen1), XYwidth)
; Now use as Q-body coil the volume that has been already reconstructed, in order to enhance SNR and CNR
; SMOOTH FILTERING
; Separate the real part of the reconstructed volume from the imaginary part, filter them with a smooth filter and compute the module
;ImagVolume=fltarr(outputRes(0),outputRes(1),outputRes(3))
;RealVolume=fltarr(outputRes(0),outputRes(1),outputRes(3))
;ImagVolume=imaginary(volumen1)
;RealVolume=real_part(volumen1)
;kernel=intarr(3)
;kernel=[5,5,5]
;ImagVolume=SmoothFilter(temporary(ImagVolume),kernel)
;RealVolume=SmoothFilter(temporary(RealVolume),kernel)
;QBodyCoil=sqrt(ImagVolume^2+RealVolume^2)
; SECOND RECONSTRUCTION
SVDthreshold=1e-2
volumen2=final_reconstruction(SVDthreshold, outputRes, data, SENSEFactor, SenseMaps, QBodyCoil)
volumen2=temporary(volumen2)*abs(ExportMask)
stop
end