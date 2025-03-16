# TODO:
# -then the function that works for vectors should be able to work for matrices with minimal correction rather than using the base version.
# change the names so that they are more descriptive
# - continue with the naming convention so that a variable can be tracked across different files
# DONE:
# -make all the variables that are row major to be column major. E.g. eachrow
"""
    This function calculates the anomalies for the temperature array and the clim/threshold.
- anom: the anomaly
- thsd: the threshold difference
- anbx (index): for the rate of onset
- anfx (index): for the rate of decline

- anbf (value): for the rate of onset
- anft (value): for the rate of decline
"""

function _anoms(sst, clim, thsh, lyd, mst, mse)
    anom = sst[mst:mse] - clim[lyd[mst:mse]]
    thsd = thsh[lyd[mst:mse]] - clim[lyd[mst:mse]]
    lnxt = lastindex(sst)
    anbx, anfx = max(1, mst - 1), min(lnxt, mse + 1)
    anbf = sst[anbx] - clim[lyd[anbx]]
    anft = sst[anfx] - clim[lyd[anfx]]
    return (; anom, anbf, anft, lnxt, thsd)
end

_categorys(anoms) = min(4, maximum(fld.(anoms.anom, anoms.thsd)))

function ronset(anoms, mst, anomfn, argfn)
    anom, anbf = anoms.anom, anoms.anbf
    nmx, fan, ngx = anomfn(anom), first(anom), argfn(anom)
    lnmx, snmx = nmx - mean((fan, anbf)), nmx - fan
    ron = mst > 1 ? /(lnmx, (ngx + 0.5)) : /(snmx, ngx)
    return ron
end

function rdecline(anoms, mse, anomfn, argfn)
    anft, lnx, nmx, lan, ngx, lnt = (anoms.anft, anoms.lnxt, anomfn(anoms.anom), last(anoms.anom), argfn(anoms.anom), length(anoms.anom))
    lnmx, snmx = (nmx - mean((lan, anft)), nmx - lan)
    rdc = mse < lnx ? /(lnmx, (lnt - ngx + 0.5)) : ngx == lnx ? /(snmx, 1.0) : /(snmx, (lnt - ngx))
    return rdc
end

# TODO: Rather than number indices, use the tuple name so we know what variable is going where. That's the reason behind NameTuples in the first place.

function matevents(msms::MarineHW, msst::MarineHeatWave, mexc::BitMatrix, mstartsxs, mendsxs)
    mni, cmi, mxi, drt, cats, rons, rdcs, vri, sdt, edt, cls, rws = ntuple(_ -> Vector(undef, length(mstartsxs)), 12)
    fullyears = unique(year.(msst.dates))
    mask = msst.mask
    x, y = size(msms.temp)
    climout, threshout = ntuple(_ -> similar(msms.temp, x, y, 366), 2)
    msms.exceed[mask, :] = mexc

    for i in eachindex(mstartsxs)
        evmeta = _eventmetrics(eachcol(msst.temp)[i], eachcol(msst.clima)[i], eachcol(msst.thresh)[i], msst.lyday, msst.anomfn, msst.argfn, msst.dates, mstartsxs[i], mendsxs[i])
        meanmets, evmets = meanmetrics(evmeta[1], evmeta[2], evmeta[3], fullyears, msst.anomfn)
        cats[i] = evmeta[4]
        anmets = annualmetrics(evmeta[1], evmeta[2], evmeta[3], evmeta[5], evmeta[6], fullyears, msst.anomfn)
        msms.temp[mask[i], :] = evmeta[7]
        msms.category[mask[i], :] = evmeta[8]
        climout[mask[i], :] = eachcol(msst.clima)[i]
        threshout[mask[i], :] = eachcol(msst.thresh)[i]

        rons[i], rdcs[i], sdt[i], edt[i] = evmeta[2], evmeta[3], evmeta[5], evmeta[6]
        rws[i] = repeat([mask[i][1]], length(evmeta[2]))
        cls[i] = repeat([mask[i][2]], length(evmeta[2]))
        mni[i], cmi[i], mxi[i], drt[i], vri[i] = evmets
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

function mhctemp(sst, sstdate::StepRange{Date,Day}, mdate::StepRange{Date,Day}, cdate::StepRange{Date,Day}; threshold=threshold, pwidth=31, win_width=5)
    sstdate, mdate, cdate = testdates(sstdate, mdate, cdate)
    testarrays(sst, sstdate)
    mhsst, mask, mlyd = subtemp(sst, sstdate, mdate)
    clsst, mask, clyd = subtemp(sst, sstdate, cdate)
    dvec = daterange(clyd, win_width)
    clima, climq = clthr(clsst, dvec, threshold, pwidth)
    return mhsst, mdate, mlyd, mask, clima, climq
end


function edetect(sst::Vector, sstdate, mdate, cdate; threshold=0.9)
    excdfn = ≥
    msst = MHTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    msms = vecarr(sst, mdate)
    mexc = exceed(msst)
    msxs, mexs = _mylabeling(mexc, min_dur, max_gap)
    evanom, rons, rdcs, cats, stdate, endate, mhwout, catout = eventmetrics(msst, msxs, mexs)
    fullyears = unique(year.(mdate))
    meanmets, evmets = meanmetrics(evanom, rons, rdcs, fullyears, msst.anomfn)
    anmets = annualmetrics(evanom, rons, rdcs, stdate, endate, fullyears, msst.anomfn)
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


function edetect(sst::Array, sstdate, mdate, cdate; threshold=0.9)
    excdfn = >=
    msst = MHTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    msms = vecarr(sst, mdate)
    mexc = exceed(msst)
    msxs, mexs = _mylabeling(mexc, min_dur, max_gap)
    msms, edf, cl, th = matevents(msms, msst, mexc, msxs, mexs)
    return msms, edf, cl, th
end

function edetect(msst::MarineHeatWave, msms::MarineHW)
    mexc = exceed(msst)
    msxs, mexs = _mylabeling(mexc, min_dur, max_gap)
    msms, edf, cl, th = matevents(msms, msst, mexc, msxs, mexs)
    return msms, edf, cl, th
end

function matevents(msms::MarineHW, msst::MarineHeatWave, mexc::BitVector, mstartsxs, mendsxs)
    evanom, rons, rdcs, cats, stdate, endate, mhwout, catout = eventmetrics(msst, mstartsxs, mendsxs)
    fullyears = unique(year.(msst.dates))
    meanmets, evmets = meanmetrics(evanom, rons, rdcs, fullyears, msst.anomfn)
    anmets = annualmetrics(evanom, rons, rdcs, stdate, endate, fullyears, msst.anomfn)
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
    msst = MHTemp(mhctemp(sst, sdate, mdate, cdate; threshold)..., excdfn, argfn, anomfn)
    msms = vecarr(sst, mdate)
    msms, edf, cl, th = edetect(msst, msms)
    return msms, edf, cl, th
end

function MarineCS(sst::Array, sdate, mdate, cdate; threshold=0.1)
    excdfn, argfn, anomfn = ≤, argmin, minimum
    msst = MCTemp(mhctemp(sst, sdate, mdate, cdate; threshold)..., excdfn, argfn, anomfn)
    msms = vecarr(sst, mdate)
    msms, edf, cl, th = edetect(msst, msms)
    return msms, edf, cl, th
end
