Murthy lab version of the Flysongsegmenter ([Arthur et al., 2013][1]).
The [original code][2] was modified by Pip Coen ([Coen et al., 2015][3]).
Now includes a [pulse type classifier][7].

See `demo_singleMicChamber.m` and `demo_multiMicChamber.m` for usage and comparison with hand-annotated song data. `demo_manualAnnotation.m` shows how to start the [manual segmenter][6].

The song for the single-mic demo is a recording from the data supplement of [Stern (2014)][4]. The hand-annotation for that song comes from the supplemental material of [Kyriacou et al. (2017)][5]. The multi-mic demo uses an unpublished recording of _Drosophila melanogaster_ and was hand-annotated using the [manual segmenter][6].

[1]: https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-11-11 "Multi-channel acoustic recording and automated analysis of Drosophila courtship songs"
[2]: https://github.com/FlyCourtship/FlySongSegmenter "fly song segmenter"
[3]: https://www.nature.com/nature/journal/v507/n7491/full/nature13131.html "Dynamic sensory cues shape song structure in _Drosophila_"
[4]: https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-12-38 "Reported Drosophila courtship song rhythms are artifacts of data analysis"
[5]: http://www.pnas.org/content/114/8/1970.abstract "Failure to reproduce period-dependent song cycles in Drosophila is due to poor automated pulse-detection and low-intensity courtship"
[6]: demo_manualAnnotation.m "demo_manualAnnotation.m"
[7]: https://github.com/postpop/pulseTypeClassifier "pulse type classifier"
