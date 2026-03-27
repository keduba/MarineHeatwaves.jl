## Four simple steps

### to a marine heatwave calculation.

Maybe even less, depending on what you want to achieve.

### 1. Create your marine extreme object.

Assuming you already loaded your temperature data (sst), you would also need your dates - sst date, climatology date, and marine heatwave (or cold spell) date.

So three of them: sst\_date, clim\_date and mhw/mcs\_date.

`mex, cixs = mextreme(sst, sstdate, mhwdate, climdate)`

Here `mex` is your marine heatwave or coldspell object which you would use for calculating the event itself in the next steps.

`cixs` contains the indices where you had non-NaN data ahead of an internal transformation. It is passed later as an argument to return a full array.

----------

### 2. Find the days where there was an event

You probably know how long in days you want as your cutoff for an event. Usually it's 5 days (but you could experiment with longer `;)`). This should be passed as `minimum_duration`. How long should be between events for it to qualify as one event? The default gap is 2 days but again, you could try reasonable values to see.

Alright, minimum\_duration = 5, maximum\_gap = 2. These are the defaults.

`sts_ends = mhlabels(mex, min_dur, max_gap)`

If your original temperature data (sst) was an array (3D), then `sts_ends` will contain a tuple of the starts, ends and pixels where events were found. It can happen that certain pixels never have events especially in cold spells.

----------

### 3. Calculate the characteristics of the events

Easy as saying

`mevts = anomsa(mex, sts_ends)`

With `mevts`, you have all the events and their characteristics which you can view directly as a table or matrix or dataframe as you prefer.

To view it immediately, as a matrix:

`mymetric(mevts)` if your original input was a single point (vector).

`mymetric(mevts, indices, sts_ends, mhwdate)`

This returns you a full matrix of the events with the pixels in which they occurred, and the event dates.

----------

### 4. And what about the metrics per pixel?

#### The mean metrics?

`mean_mets = meanmetrics(mevt, mhwdate)`

#### Or the annual metrics?

`ann_mets = annualmetrics(mevt, mhwdate, sts_ends)`

#### And the trends over the given years?

`trmets = trend(ann_mets)`

You can also obtain these values as one object with:

`all_mets = mhmetrics(mevt, mhwdate, sts_ends)`

----------

## View your results

`mymetric(mevt)`

`mymetric(all_mets)`
