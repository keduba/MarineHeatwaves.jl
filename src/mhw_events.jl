

function matevents(msms::MarineHW, msst::MarineHeatWave, mexc::BitMatrix, mstartsxs, mendsxs)
    mni, cmi, mxi, drt, cats, rons, rdcs, vri, sdt, edt, cls, rws = ntuple(_ -> Vector(undef, length(mstartsxs)), 12)
    fullyears = unique(year.(msst.dates))
    mask = msst.mask
    x, y = size(msms.temp)
    climout, threshout = ntuple(_ -> similar(msms.temp, x, y, 366), 2)
    msms.exceed[mask, :] = mexc
    for i in eachindex(mstartsxs)
        evmeta = _eventmetrics(eachrow(msst.temp)[i], eachrow(msst.clima)[i], eachrow(msst.thresh)[i], msst.lyday, msst.anomfn, msst.argfn, msst.dates, mstartsxs[i], mendsxs[i])
        cats[i] = evmeta[4]
        msms.temp[mask[i], :] = evmeta[7]
        msms.category[mask[i], :] = evmeta[8]
        climout[mask[i], :] = eachrow(msst.clima)[i]
        threshout[mask[i], :] = eachrow(msst.thresh)[i]
        meanmets, evmets = meanmetrics(evmeta[1], evmeta[2], evmeta[3], fullyears)
        rons[i], rdcs[i], sdt[i], edt[i] = evmeta[2], evmeta[3], evmeta[5], evmeta[6]
        rws[i] = repeat([mask[i][1]], length(evmeta[2]))
        cls[i] = repeat([mask[i][2]], length(evmeta[2]))
        mni[i], cmi[i], mxi[i], drt[i], vri[i] = evmets
        anmets = annualmetrics(evmeta[1], evmeta[2], evmeta[3], evmeta[5], evmeta[6], fullyears)
        trmets = trends(anmets)
        for m in keys(msms.annuals)
            msms.annuals[m][mask[i], :] = anmets[m]
            msms.means[m][mask[i]] = meanmets[m]
            msms.trends[m][mask[i]] = trmets.trends[m]
            msms.pvalues[m][mask[i]] = trmets.pvalues[m]
            msms.pmetrics[m][mask[i]] = trmets.pmetrics[m]
        end
    end

    edf = (MeanInt=reduce(vcat, mni), CumInt=reduce(vcat, cmi), MaxInt=reduce(vcat, mxi), Duration=reduce(vcat, drt), Category=reduce(vcat, cats), ROnset=reduce(vcat, rons), RDecline=reduce(vcat, rdcs), VarInt=reduce(vcat, vri), StartDate=reduce(vcat, sdt), EndDate=reduce(vcat, edt), Xind=reduce(vcat, rws), Yind=reduce(vcat, cls))

    return msms, edf, climout, threshout
end

function mhctemp(sst, sstdate::StepRange{Date,Day}, mdate::StepRange{Date,Day}, cdate::StepRange{Date,Day}; threshold=threshold)
    sstdate, mdate, cdate = testdates(sstdate, mdate, cdate)
    testarrays(sst, sstdate)
    mhsst, mask, mlyd = subtemp(sst, sstdate, mdate)
    clsst, mask, clyd = subtemp(sst, sstdate, cdate)
    dvec = daterange(clyd, winwidth)
    clima, climq = clthr(clsst, dvec, threshold)
    return mhsst, mdate, mlyd, mask, clima, climq
end


function edetect(sst::Vector, sstdate, mdate, cdate; threshold=0.9)
    excdfn = ≥
    mhwin = MHTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    msms = vecarr(sst, mdate)
    mexc = exceed(mhwin)
    msxs, mexs = _mylabeling(mexc)
    evanom, rons, rdcs, cats, stdate, endate, mhwout, catout = eventmetrics(mhwin, msxs, mexs)
    fullyears = unique(year.(mdate))
    meanmets, evmets = meanmetrics(evanom, rons, rdcs, fullyears)
    anmets = annualmetrics(evanom, rons, rdcs, stdate, endate, fullyears)
    trmets = trends(anmets)
    msms.exceed .= mexc
    msms.temp .= mhwout
    msms.category .= catout
    # Annuals
    for m in keys(msms.annuals)
        msms.annuals[m][begin:end] = anmets[m]
        msms.means[m] .= meanmets[m]
        msms.trends[m][1] = trmets.trends[m]
        msms.pvalues[m][1] = trmets.pvalues[m]
        msms.pmetrics[m][1] = trmets.pmetrics[m]
    end
    edf = (MeanInt=evmets[1], CumInt=evmets[2], MaxInt=evmets[3], Duration=evmets[4], Category=cats, ROnset=rons, RDecline=rdcs, VarInt=evmets[5], StartDate=stdate, EndDate=endate)
    return msms, edf, mhwin.clima, mhwin.thresh
end


function edetect(sst::Array, sstdate, mdate, cdate; threshold=0.9)
    excdfn = >=
    mhwin = MHTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    msms = vecarr(sst, mdate)
    mexc = exceed(mhwin)
    msxs, mexs = _mylabeling(mexc)
    msms, edf, cl, th = matevents(msms, mhwin, mexc, msxs, mexs)
    return msms, edf, cl, th
end

function edetect(mhwin::MarineHeatWave, msms::MarineHW)
    mexc = exceed(mhwin)
    msxs, mexs = _mylabeling(mexc)
    msms, edf, cl, th = matevents(msms, mhwin, mexc, msxs, mexs)
    return msms, edf, cl, th
end

function matevents(msms::MarineHW, msst::MarineHeatWave, mexc::BitVector, mstartsxs, mendsxs)
    evanom, rons, rdcs, cats, stdate, endate, mhwout, catout = eventmetrics(msst, mstartsxs, mendsxs)
    fullyears = unique(year.(msst.dates))
    meanmets, evmets = meanmetrics(evanom, rons, rdcs, fullyears)
    anmets = annualmetrics(evanom, rons, rdcs, stdate, endate, fullyears)
    trmets = trends(anmets)
    msms.exceed .= mexc
    msms.temp .= mhwout
    msms.category .= catout
    # Annuals
    for m in keys(msms.annuals)
        msms.annuals[m][begin:end] = anmets[m]
        msms.means[m] .= meanmets[m]
        msms.trends[m][1] = trmets.trends[m]
        msms.pvalues[m][1] = trmets.pvalues[m]
        msms.pmetrics[m][1] = trmets.pmetrics[m]
    end
    edf = (MeanInt=evmets[1], CumInt=evmets[2], MaxInt=evmets[3], Duration=evmets[4], Category=cats, ROnset=rons, RDecline=rdcs, VarInt=evmets[5], StartDate=stdate, EndDate=endate)
    return msms, edf, msst.clima, msst.thresh
end

function MarineHW(sst::Array, sdate, mdate, cdate; threshold=0.9)
    excdfn, argfn, anomfn = ≥, argmax, maximum
    mhwin = MHTemp(mhctemp(sst, sdate, mdate, cdate; threshold)..., excdfn, argfn, anomfn)
    msms = vecarr(sst, mdate)
    msms, edf, cl, th = edetect(mhwin, msms)
    return msms, edf, cl, th
end

function MarineCS(sst::Array, sdate, mdate, cdate; threshold=0.1)
    excdfn, argfn, anomfn = ≤, argmin, minimum
    mhwin = MCTemp(mhctemp(sst, sdate, mdate, cdate; threshold)..., excdfn, argfn, anomfn)
    msms = vecarr(sst, mdate)
    msms, edf, cl, th = edetect(mhwin, msms)
    return msms, edf, cl, th
end
