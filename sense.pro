pro HomodyneRecon, data,outputRes=outputRes, filter=filter
print, "Homodyne filtering"
  filter=float(filter)
  filter=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.)
  mask=filter;float(erode(filter, replicate(1.,10,10)))
  a=shift(reverse(shift(reverse(filter,1),1,0),2),0,1)
  b=filter+a
  ind=where(b eq 2, count)
    a=fltarr(outputRes(0),outputRes(1))
    a(ind)=1
    a=a+shift(a,10,0)
    a=a+shift(a,0,10)
    ind=where(a gt 0)
  filter(ind)=0.7
  filter=filter*mask
  filter=convol(filter, FilterGen(30,19),/EDGE_t)

     weight=fltarr(1,outputRes(1))
     weight(0,*)=exp(-(findgen(outputRes(1))-outputRes(1)/2.)^2./50000)
     weight=rebin(weight, outputRes(0), outputRes(1))

     a=fltarr(outputRes(0),1)
     a(*,0)=exp(-(findgen(outputRes(0))-outputRes(0)/2.)^2./50000)
     a=rebin(a, outputRes(0), outputRes(1))

     weight *= a

  filter=filter*weight
  tvscl, filter
  filter=shift(filter,outputRes(0)/2.,outputRes(1)/2.)
  filter=rebin(filter, outputRes(0),outputRes(1),outputRes(3))

  data=fft(fft(data, 1,dimension=1), 1, dimension=2); pasa al espacio-K
  for i=0, outputRes(4)-1 do begin
    inter=data(*,*,*,i)*filter
    data(*,*,*,i)=inter
  endfor
  data=fft(fft(data, -1,dimension=1), -1, dimension=2); pasa al espacio real
end

Function Generate_KSpace_mask, maskKspace
    ;maskKspace=(total(total(maskKspace,2),2),2)
    dim=size(maskKspace, /dimension)
    KSpace_mask=shift(maskKspace, -dim(0)/2.,-dim(1)/2.)
    KSpace_mask=smooth(KSpace_mask,5)
    ind=where(KSpace_mask gt 0)
    KSpace_mask(ind)=1.
    KSpace_mask=erode(KSpace_mask,replicate(1.,5,5))
    KSpace_mask=shift(KSpace_mask, dim(0)/2.,dim(1)/2.)
    return, KSpace_mask
end

pro SENSE
  argument=command_line_args(count=count)

  dir=argument(0)
;
dir='C:\BrisaPatients\BRISA085\'
;  stop
loadct,0
  if file_test(dir) then begin
     genVideo=1
     FOVReduction=0
     Recursive=0
     mask=1
     fullSamplingMask=1
     WriteData=1
     RingingFilterOrder=50
     XYshift=[10,-5];X=Negative Left must be multiplos de 6; Y=Negative Anterior must be multiplos de 4

     if file_test(dir+'ReconConfig.txt') then begin
          params=((read_ascii(dir+'ReconConfig.txt')).field1)(*,1)
          genVideo=params(0)
          FOVReduction=params(1)
          Recursive=params(2)
          mask=params(3)
          fullSamplingMask=params(4)
          XYshift=params(5:6)
     endif

     ReadingTime=0.0
     ReconTime=0.0
     SensesitivitiesTime=0.0

     ReadingTime=SYSTIME(/SECONDS )
     RefProto=ReadProtocol(dir+'Ref_XL_Torso.txt');dialog_pickfile(filter='*.txt'))
     ScanProto=ReadProtocol(dir+'CORONARIAS+EDEMA.txt')

     ;Fractional SENSE correction
     RealImageFOV=[ScanProto.FOV(1),ScanProto.FOV(2),ScanProto.FOV(0)]

     ScanProto.FOV(1:2)=ScanProto.FOV(1:2)*ceil(ScanProto.SENSE)/ScanProto.SENSE

     TotalSlices  = round(ScanProto.FOV(0)/ScanProto.VOXELSIZE(0))
     SlicesCenter = 0 ; positive is food direction negative head direction

     RR_interval=((60./float(ScanProto.CardFrequency))*1000.)/float(ScanProto.CardPhases)

     CoilHeaderCPX=file_search(dir+'cpx*.list')
     CoilDataCPX=file_search(dir+'cpx*.data')

     DataHeader=file_search(dir+'raw*.list')
     DataData=file_search(dir+'raw*.data')

     Coils=read_cpx(CoilHeaderCPX, CoilDataCPX, QBodyCoil, NCoils)
     outputRes=float([640,640, NCoils, 20, 1])
     ; FOV shift due to different FOV sampling scheme for fractional SENSE factors
     FOVShift=(ceil(ScanProto.SENSE)-ScanProto.SENSE)/ceil(ScanProto.SENSE)
     ScanProto.SENSE=ceil(ScanProto.SENSE)
     ReadingTime=SYSTIME(/SECONDS )-ReadingTime

     SensesitivitiesTime=SYSTIME(/SECONDS )
     FOVCorrection=scanproto.FOV
     data=read_raw_list_by_TFEFactor(DataHeader, DataData, outputRes=outputRes, slice=SlicesCenter,ScanProto=ScanProto, noise,SENSEFactor=ScanProto.SENSE, mask=maskKspace,TimeMask=TimeMask);RR_interval=RR_interval,
     ;prepare data for reconstruction
     FOVCorrection=scanproto.FOV/FOVCorrection

     window,1, xsize=outputRes(0)*4., ysize=outputRes(1)*4., title=dir
     for i=0, outputRes(4)-1 do tvscl,shift(alog(abs(data(*,*,0,outputRes(3)/2,i))),outputRes(0)/2,outputRes(1)/2) ,i

     if total(ScanProto.SENSE) lt 2.0 then ScanProto.SENSE= [3.0,2.0]

;NoiseDeco,SenseMaps,data, noise, outputRes=outputRes

;    dummy=NormData(data, maskKspace, outputRes=outputRes,TimeMask=TimeMask)

     maskKspace=Generate_KSpace_mask(maskKspace)
     vectorX=(where(maskKspace(*,0) eq 0))
     vectorY=(where(maskKspace(0,*) eq 0))
     radio=float([vectorX(0),vectorY(0)])+1
     LowerRadio=radio-[3.,3.]
     HalfScan=[vectorX(n_elements(vectorX)-1),vectorY(n_elements(vectorY)-1)]
     LowerHalfScan=HalfScan+[3.,3.]

     maskKspace=fltarr(outputRes(0),outputRes(1))
     filter=fltarr(outputRes(0),outputRes(1))
     for j=0, outputRes(1)-1 do begin
       for i=0, outputRes(0)-1 do begin
            radius=((i-outputRes(0)/2.)/radio(0))^2.+((j-outputRes(1)/2.)/radio(1))^2.
            if radius le 1 then maskKspace(i,j)=1.0
            LowerRadius=((i-outputRes(0)/2.)/LowerRadio(0))^2.+((j-outputRes(1)/2.)/LowerRadio(1))^2.
            if LowerRadius le 1 then filter(i,j)=1.0
       endfor
     endfor

     maskKspace=shift(maskKspace,outputRes(0)/2.,outputRes(1)/2.)
     maskKspace(radio(0):HalfScan(0),*)=0
     maskKspace(*,radio(1):HalfScan(1))=0

     filter=shift(filter,outputRes(0)/2.,outputRes(1)/2.)
     filter(LowerRadio(0):LowerHalfScan(0),*)=0
     filter(*,LowerRadio(1):LowerHalfScan(1))=0

     HomodineFilter=maskKspace

     filter=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.)
     filter=(convol(filter, FilterGen(25,19),/EDGE_t))

     maskKspace=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.)

     weight=fltarr(1,outputRes(1))
     weight(0,*)=exp(-(findgen(outputRes(1))-outputRes(1)/2.)^2./40000)
     weight=rebin(weight, outputRes(0), outputRes(1))

     a=fltarr(outputRes(0),1)
     a(*,0)=exp(-(findgen(outputRes(0))-outputRes(0)/2.)^2./40000)
     a=rebin(a, outputRes(0), outputRes(1))

     weight *= a

     filter=filter*shift(HomodineFilter, -outputRes(0)/2.,-outputRes(1)/2.)*weight
     weight=0

     filter=reform(shift(filter, -outputRes(0)/2.,-outputRes(1)/2.),  outputRes(0),outputRes(1),1,1)
     filter=rebin(filter,outputRes(0),outputRes(1),1,outputRes(3))

     for i=0, outputRes(2)-1 do begin
        for j=0, outputRes(4)-1 do begin
          data(*,*,i,*,j)=data(*,*,i,*,j)*filter
        endfor
     endfor
     filter=reform(filter)

     mask=fltarr(outputRes(0),outputRes(1),outputRes(2),1,outputRes(4))
     ind=where(abs(data(*,*,*,0,*)), count)
     mask(ind)=1.

     mask=total(mask,5)
     ind=where(abs(mask), count)

     StaticKSpaceData=total(data,5)
     for i=0, outputRes(3)-1 do begin
        inter=Reform(StaticKSpaceData(*,*,*,i))
        inter(ind)=inter(ind)/mask(ind)
        StaticKSpaceData(*,*,*,i)=reform(inter,outputRes(0),outputRes(1),outputRes(2),1)
     end

     StaticOriginalKSpaceDataIndex=where(abs(reform(StaticKSpaceData(*,*,0,*))) gt 0, countOriginal)

     Data=Data(ScanProto.SENSE(0)*findgen(outputRes(0)/(ScanProto.SENSE(0))),*,*,*,*)
     Data=Data(*,ScanProto.SENSE(1)*findgen(outputRes(1)/(ScanProto.SENSE(1))),*,*,*)
     for i=0, outputRes(4)-1 do tvscl, shift(alog(abs(data(*,*,0,outputRes(3)/2,i))),-outputRes(0)/ScanProto.SENSE(0)/4. ,-outputRes(1)/ScanProto.SENSE(1)/4. ),i
     imagen=tvrd()
     WRITE_PNG, dir+'samplingAlongPhases.png', Imagen


     ChannelSNR=total(total(abs(data(0,0,*,*,*)),4),4)/outputRes(3)/outputRes(4)
     ind=reverse(sort(ChannelSNR))
     ChannelSelection=where(ChannelSNR gt mean(ChannelSNR(ind(0:2.*outputRes(2)/3.)))-2.0*stddev(ChannelSNR(ind(0:outputRes(2)/2.))), NumSelectedChannels)

     if NumSelectedChannels lt 8 then NumSelectedChannels=outputRes(2)

     if NumSelectedChannels lt outputRes(2) then begin
          print, ChannelSelection
          ChannelSNR=0
          maskKspace=maskKspace(*,*,ChannelSelection,*,*)
          data=data(*,*,ChannelSelection,*,*)
          StaticKSpaceData=StaticKSpaceData(*,*,ChannelSelection,*)
     endif

     fullSamplingMask=1
     SensesitivitiesTime=SYSTIME(/SECONDS )
     SenseMaps=SenseMapsGen(Coils, QBodyCoil, RefProto, ScanProto, outputRes=outputRes, QBThreshold=5e-2, mask=fullSamplingMask,FOVCorrection=FOVCorrection)
     SensesitivitiesTime=SYSTIME(/SECONDS )-SensesitivitiesTime

     if (NumSelectedChannels lt outputRes(2)) then begin
        SenseMaps=temporary(SenseMaps(*,*,*,ChannelSelection))
        outputRes(2)=NumSelectedChannels
     end

     Coils=0
     QBodyCoil=abs(QBodyCoil)

     ReconTime=SYSTIME(/SECONDS)

     Data=SENSE_all_in_one(Data, SenseMaps,outputRes=outputRes, SENSEFactor=ScanProto.SENSE,SVDthreshold=1e-3,QBodyCoil=QBodyCoil,fullSamplingMask=fullSamplingMask)
     ReconTime=SYSTIME(/SECONDS)-ReconTime


;     HomodyneRecon, data,outputRes=outputRes, filter=HomodineFilter

     filename=dir+'DICOM\IM_0001'
     info=read_dicom_info(filename)
     help, info,/struct
     if keyword_set(WriteData) then begin
            print, 'Resizing to the final image resolution'
            FinalImageFOV=[ScanProto.FOV(1), Scanproto.FOV(2), Scanproto.FOV(0)]
            NewDim=round(FinalImageFOV/[info.pixelsize(0),info.pixelsize(2),info.pixelsize(1)])
            NewDim=[NewDim,fix(outputRes(4))]
            Inter=fltarr(NewDim(0),NewDim(1),NewDim(2),outputRes(4))
            SlicePos=[round((NewDim(0)-outputRes(0))/2.),round((NewDim(1)-outputRes(1))/2.),round((NewDim(2)-outputRes(3))/2.)]

            for i=0, outputRes(4)-1 do begin
                a=complexarr(NewDim(0),NewDim(1),NewDim(2))
                b=fft(reform(data(*,*,*,i)),-1)
                b=shift(b,-round(outputRes(0)/2.),-round(outputRes(1)/2.),-round(outputRes(3)/2.))
                a(SlicePos(0):SlicePos(0)+outputRes(0)-1,SlicePos(1):SlicePos(1)+outputRes(1)-1,$
                  SlicePos(2):SlicePos(2)+outputRes(3)-1)=b
                ;a=a*Filter
                a=shift(a, round(NewDim(0)/2.),round(NewDim(1)/2.),round(NewDim(2)/2.))
                a=fft(a,1)
                Inter(*,*,*,i)=abs(a)
                a=0
                b=0
            endfor

            data=fltarr(info.npixels(0), info.npixels(2), info.npixels(1),outputRes(4))

            if NewDim(0) lt info.npixels(0) then begin
             SlicePos=[round((info.npixels(0)-NewDim(0))/2.),round((NewDim(1)-info.npixels(2))/2.),round((info.npixels(1)-NewDim(2))/2.)]
               data(SlicePos(0):SlicePos(0)+NewDim(0)-1,*,SlicePos(2):SlicePos(2)+NewDim(2)-1,*)=$
                    Inter(*,SlicePos(1):SlicePos(1)+info.npixels(2)-1,*,*)
            endif else begin
              SlicePos=[round((NewDim(0)-info.npixels(0))/2.),round((NewDim(1)-info.npixels(2))/2.),round((info.npixels(1)-NewDim(2))/2.)]
              data(*,*,SlicePos(2):SlicePos(2)+NewDim(2)-1,*)=$
                    Inter(SlicePos(0):SlicePos(0)+info.npixels(0)-1,SlicePos(1):SlicePos(1)+info.npixels(2)-1,*,*)
            endelse
            Inter=0
            data=abs(data)
            maxDynImage=max(data(*,*,NewDim(2)/2-3:NewDim(2)/2+3,*))

            data=data<(maxDynImage)

            data=((data*4000./(maxDynImage))<4000.)>0

            openw,1, dir+'time.txt'
            printf,1,  "ReadingTime  SensesitivitiesTime   ReconTime"
            printf,1,  ReadingTime, SensesitivitiesTime, ReconTime
            close,1
            help, data
;            write_dicom_volume, data, dir
            View3D, dir=dir, data=data

       endif else begin
            data=abs(data)
            maxDynImage=max(data(*,*,outputRes(3)/2-3:outputRes(3)/2+3,*))

            data=data<(maxDynImage)

            data=((data*4000./maxDynImage)<4000.)>0

            data=reverse(data,1)
            data=reverse(data,2)
            View3D, dir=dir,data=data
       endelse
  endif
end