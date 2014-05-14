import numpy
import pylab
import pandas
import matplotlib
import cellh5
import cellh5_analysis
import re
import csv
from itertools import izip_longest

from PIL import Image as PImage

from estimator import HMMConstraint, HMMAgnosticEstimator, normalize

from sklearn import hmm

class IgorsClass(cellh5_analysis.CellH5Analysis):
    cmap3 = matplotlib.colors.ListedColormap(map(lambda x: cellh5_analysis.hex_to_rgb(x), 
                                                   ['#FFFFFF', 
                                                    '#11FFF11', 
                                                    '#FF0000',
                                                    '#0000FF',
                                                    '#000000']), 'cmap3')
    
    cmap4 = matplotlib.colors.ListedColormap(map(lambda x: cellh5_analysis.hex_to_rgb(x), 
                                                   ['#FFFFFF', 
                                                    '#11FFF11', 
                                                    '#FF0000',
                                                    '#0000FF',
                                                    '#FFFF00',
                                                    '#000000']), 'cmap4')
    cmap5 = matplotlib.colors.ListedColormap(map(lambda x: cellh5_analysis.hex_to_rgb(x), 
                                                   ['#FFFFFF', 
                                                    '#11FFF11', 
                                                    '#FF0000',
                                                    '#0000FF',
                                                    '#05AC04',
                                                    '#FFFF00',
                                                    '#000000']), 'cmap5')

    def setup_hmm(self, hmm_n_classes, hmm_n_obs, hmm_constraint_file):
        self.hmm_n_classes = hmm_n_classes
        self.hmm_n_obs = hmm_n_obs
        self.hmm_constraint_file = hmm_constraint_file
        
        constraints = HMMConstraint(self.hmm_constraint_file)
        
        transmat = numpy.array([
                                [1.0, 0.1, 0.0, 0.0, 0.0,],
                                [0.0, 1.0, 0.1, 0.0, 0.0,],
                                [0.0, 0.0, 1.0, 0.1, 0.0,],
                                [0.0, 0.0, 0.0, 10 , 0.1,],
                                [0.1, 0.0, 0.0, 0.0, 10 ,],
                                ])
        transmat = normalize(transmat, axis=1, eps=0.00001 )
        
        assert transmat.shape[0] == self.hmm_n_classes
        assert transmat.shape[1] == self.hmm_n_classes
        
        est = HMMAgnosticEstimator(self.hmm_n_classes, 
                                   transmat, 
                                   numpy.ones((self.hmm_n_classes, self.hmm_n_obs)), 
                                   numpy.ones((self.hmm_n_classes, )))
        
        est.constrain(constraints)
        self.hmm = hmm.MultinomialHMM(n_components=est.nstates, transmat=transmat, startprob=est.startprob, init_params="")
        self.hmm._set_emissionprob(est.emis)
        
    def get_mitotic_timing(self):
        result = {}
        pattern = r".*1+(?P<mito>(2+3+))4+.*"
        for _, (plate_name, w, p, t1, t2, track_ids, track_labels) in self.mapping[['Plate', 
                                                                                    'Well', 
                                                                                    'Site', 
                                                                                    'Gene Symbol',
                                                                                    'siRNA ID',
                                                                                    'Event track ids',
                                                                                    'Event HMM track labels']].iterrows(): 
            timings = []
            for _, track_label in zip(track_ids, track_labels):
                track_str = "".join(map(str, list(track_label))) 
                mito_re = re.search(pattern, track_str)
                if mito_re is not None:
                    timings.append((mito_re.end("mito") - mito_re.start("mito")) * self.time_lapse[plate_name])
                    
            result[(w,p,t1,t2)] = timings
        return result

            

        
    def select_tracks(self, output_name):
        all_selected_tracks = []
        
        pattern = r".*1+(?P<mito_2_mito>(2+3+(4|5|1){100,}))2+3+.*"
        result_timing = {}
        for _, (plate_name, w, p, t1, t2, track_ids, track_labels) in self.mapping[['Plate', 
                                                                                    'Well', 
                                                                                    'Site', 
                                                                                    'Gene Symbol',
                                                                                    'siRNA ID',
                                                                                    'Event track ids',
                                                                                    'Event HMM track labels']].iterrows(): 
            selected_track = []
            timings = []
            for track_id, track_label in zip(track_ids, track_labels):
                track_str = "".join(map(str, list(track_label))) 
                mito_re = re.search(pattern, track_str)
                
                if mito_re is not None:
                    print "IN ", track_str
                    selected_track.append(track_label)
                    print mito_re.start("mito_2_mito"), "-->", mito_re.end("mito_2_mito")
                    timings.append((mito_re.end("mito_2_mito") - mito_re.start("mito_2_mito")) * self.time_lapse[plate_name]/60)
                    
            all_selected_tracks.append(selected_track)
            result_timing[(w, p, t1, t2)] = timings
        self.mapping[output_name] = pandas.Series(all_selected_tracks)
        return result_timing
      
        
    
    def select_tracks_G1(self, output_name):
        all_selected_tracks_G1 = []
        
        pattern = r".*1+2+3+(?P<mito_unknown>(4+))$"
        result_timing_G1 = {}
        for _, (plate_name, w, p, t1, t2, track_ids, track_labels) in self.mapping[['Plate', 
                                                                                    'Well', 
                                                                                    'Site', 
                                                                                    'Gene Symbol',
                                                                                    'siRNA ID',
                                                                                    'Event track ids',
                                                                                    'Event HMM track labels']].iterrows(): 
            selected_track_G1 = []
            timings = []
            for track_id, track_label in zip(track_ids, track_labels):
                track_str = "".join(map(str, list(track_label))) 
                G1_re = re.search(pattern, track_str)
                
                if G1_re is not None:
                    print "IN ", track_str
                    selected_track_G1.append(track_label)
                    print G1_re.start("mito_unknown"), "-->", G1_re.end("mito_unknown")
                    timings.append((G1_re.end("mito_unknown") -G1_re.start("mito_unknown")) * self.time_lapse[plate_name]/60)
                    
            all_selected_tracks_G1.append(selected_track_G1)
            result_timing_G1[(w, p, t1, t2)] = timings
        self.mapping[output_name] = pandas.Series(all_selected_tracks_G1)
        return result_timing_G1
    
    def select_tracks_S(self, output_name):
        all_selected_tracks_S = []
        
        pattern = r".*4+(?P<S_phase>(5){10,})1+.*"
        result_timing_S = {}
        for _, (plate_name, w, p, t1, t2, track_ids, track_labels) in self.mapping[['Plate', 
                                                                                    'Well', 
                                                                                    'Site', 
                                                                                    'Gene Symbol',
                                                                                    'siRNA ID',
                                                                                    'Event track ids',
                                                                                    'Event HMM track labels']].iterrows(): 
            selected_track_S = []
            timings_S = []
            for track_id, track_label in zip(track_ids, track_labels):
                track_str = "".join(map(str, list(track_label))) 
                S_re = re.search(pattern, track_str)
                
                if S_re is not None:
                    print "IN ", track_str
                    selected_track_S.append(track_label)
                    print S_re.start("S_phase"), "-->", S_re.end("S_phase")
                    timings_S.append((S_re.end("S_phase") - S_re.start("S_phase")) * self.time_lapse[plate_name]/60)
                    
            all_selected_tracks_S.append(selected_track_S)
            result_timing_S[(w, p, t1, t2)] = timings_S
        self.mapping[output_name] = pandas.Series(all_selected_tracks_S)
        return result_timing_S
    
    def select_tracks_G2(self, output_name):
        all_selected_tracks_G2 = []
        
        pattern = r".*5+(?P<G2_phase>1+)2+3+.*"
        result_timing_G2 = {}
        for _, (plate_name, w, p, t1, t2, track_ids, track_labels) in self.mapping[['Plate', 
                                                                                    'Well', 
                                                                                    'Site', 
                                                                                    'Gene Symbol',
                                                                                    'siRNA ID',
                                                                                    'Event track ids',
                                                                                    'Event HMM track labels']].iterrows(): 
            selected_track_G2 = []
            timings_G2 = []
            for track_id, track_label in zip(track_ids, track_labels):
                track_str = "".join(map(str, list(track_label))) 
                G2_re = re.search(pattern, track_str)
                
                if G2_re is not None:
                    print "IN ", track_str
                    selected_track_G2.append(track_label)
                    print G2_re.start("G2_phase"), "-->", G2_re.end("G2_phase")
                    timings_G2.append((G2_re.end("G2_phase") - G2_re.start("G2_phase")) * self.time_lapse[plate_name]/60)
                    
            all_selected_tracks_G2.append(selected_track_G2)
            result_timing_G2[(w, p, t1, t2)] = timings_G2
        self.mapping[output_name] = pandas.Series(all_selected_tracks_G2)
        return result_timing_G2
    
    
    
    def read_secondary_channel_feature(self, output_name, feature_name='n2_avg'):
        all_secondary_feature = []
        
        
        for _, (plate_name, w, p, t1, t2, track_ids, track_labels) in self.mapping[['Plate', 
                                                                                    'Well', 
                                                                                    'Site', 
                                                                                    'Gene Symbol',
                                                                                    'siRNA ID', 
                                                                                    'Event track ids',
                                                                                    'Event HMM track labels']].iterrows(): 
            sec_feat = []
            ch5_file_handle = self.cellh5_handles[plate_name]
            ch5_pos = ch5_file_handle.get_position(w, str(p))
            
            feature_idx_for_mean_signal = ch5_pos.get_object_feature_idx_by_name('secondary__expanded', feature_name)
            
            for track_id, track_label in zip(track_ids, track_labels):
                PCNA_signal = ch5_pos.get_object_features('secondary__expanded')[track_id][feature_idx_for_mean_signal]
                sec_feat.append(PCNA_signal)
            all_secondary_feature.append(sec_feat)
            
        self.mapping[output_name] = pandas.Series(all_secondary_feature)               
        
    def generate_gallery_class_images(self):
        for _, (plate_name, w, p, t1, t2, event_ids, track_ids, track_labels) in self.mapping[['Plate', 
                                                                            'Well', 
                                                                            'Site', 
                                                                            'Gene Symbol',
                                                                            'siRNA ID', 
                                                                            'Event ids',
                                                                            'Event track ids',
                                                                            'Event HMM track labels']].iterrows(): 
            ch5_file_handle = self.cellh5_handles[plate_name]
            ch5_pos = ch5_file_handle.get_position(w, str(p))
            all_img_tracks = []
            print 'Generate galleries for', w, p, ":::  ",
            for ii, event in enumerate(track_labels):
                img_track = []
                for e in event:
                    img_track.append(ch5_pos.get_gallery_image(e))
                img_track = numpy.concatenate(img_track,1)
                all_img_tracks.append(img_track)
                print ii,
                
            
            print all_img_tracks[0].shape
            
            h = all_img_tracks[0].shape[0] * len(all_img_tracks)
            w_ = numpy.max(map(lambda xx: xx.shape[1], all_img_tracks))
            result_img = numpy.zeros((h, w_, 3), numpy.uint8)
            
            for img_idx, img in enumerate(all_img_tracks):
                print ' deeebug'
                print img.shape
                print result_img[img_idx*img.shape[0]:((img_idx+1)*img.shape[0]), 0:img.shape[1], :].shape
                result_img[img_idx*img.shape[0]:((img_idx+1)*img.shape[0]), 0:img.shape[1], :] = img
            
            img = PImage.fromarray(result_img)
            img.save(self.output("galleries_%s_%s_%s_%s.png" % (w,p,t1,t2)))
    
            
            
            
    
    def combine_classifiers(self, output_name):
        all_combined_classes = []
        for _, (plate_name, w, p, t1, t2, track_ids, track_labels) in self.mapping[['Plate', 
                                                                            'Well', 
                                                                            'Site', 
                                                                            'Gene Symbol',
                                                                            'siRNA ID', 
                                                                            'Event track ids',
                                                                            'Event track labels']].iterrows(): 
            combined_classes = []
            ch5_file_handle = self.cellh5_handles[plate_name]
            ch5_pos = ch5_file_handle.get_position(w, str(p))
            
            for track_id, track_label in zip(track_ids, track_labels):
                h2b_class = track_label.copy()
                pnas_class = ch5_pos.get_class_prediction('secondary__expanded')[track_id]['label_idx'] + 1
                
#                 pnas_class_tmp = pnas_class.copy()
#                 h2b_class_tmp = h2b_class.copy()
                
                inter_idx = h2b_class == 1
                pnas_class[pnas_class==3] = 1
                pnas_class[pnas_class==2]+=2
                
                combined_class = h2b_class
                combined_class[inter_idx] = pnas_class[inter_idx]
                
                combined_classes.append(combined_class)
#                 print "*"*10
#                 for h, p, c in zip(h2b_class_tmp, pnas_class_tmp, combined_class):
#                     print h, p, "->", c
            all_combined_classes.append(combined_classes)
            
        self.mapping[output_name] = pandas.Series(all_combined_classes)
    
if __name__ == "__main__":
    # Does alter the main appearance of the plots (such as background grids, etc)
    import seaborn
    
    print 'Start'
    plate_name = '140508'
    ch5_main_file = "/Volumes/groups/gerlich/members/Igor Gak/FINAL_SS_Plate1/analysis/hdf5/_all_positions.ch5"
    pos_mapping_file = "/Volumes/Seagate Backup Plus Drive/FINAL_HDF5/FINAL.txt"
    hmm_xml_constraint_file = "graph_igor_pnas_h2b_combi.xml"
    
    pm = IgorsClass('igor_test', 
                        {plate_name: pos_mapping_file}, 
                        {plate_name: ch5_main_file}, 
                        sites=(1,),
                       rows=("D","C","F",), 
                       cols=(4, 5,)
                    )
    pm.read_events(5,99999)
    pm.track_full_events()
    pm.combine_classifiers("Event labels combined")
    pm.setup_hmm(5, 4, hmm_xml_constraint_file)
    pm.predict_hmm("Event labels combined")
    
    pm.plot_track_order_map(["Event track labels", "Event labels combined", 'Event HMM track labels'], 
                             [IgorsClass.cmap3, IgorsClass.cmap4, IgorsClass.cmap5])
    
    
#    pm.read_tertiary_channel_feature('Cyclin values for HMM tracks')
#    pm.generate_gallery_class_images()
    
#    pm.event_curves('Event HMM track labels', 
#                           "secondary__expanded",
#                           "n2_avg",
#                           IgorsClass.cmap5,
#                           (0, 2000),
#                           (0, 2),
#                           0,
#                           )
    
    mito_timings = pm.get_mitotic_timing()
    
    mito2mito_timings = pm.select_tracks('Selected HMM track labels')
    
#    pm.plot_track_order_map(['Selected HMM track labels'], [IgorsClass.cmap5])
    
    mito_unknown_timings = pm.select_tracks_G1('Selected HMM track labels')
    
    S_phase_timings = pm.select_tracks_S('Selected HMM track labels')
    
    G2_phase_timings = pm.select_tracks_G2('Selected HMM track labels')
    
#    pm.plot_track_order_map(['Selected HMM track labels'], [IgorsClass.cmap5])
    
#    def median(S_phase_timings.values()):

#    srtd = sorted(S_phase_timings.values()) 
#    mid = len(S_phase_timings.values())/2   

#    if len(S_phase_timings.values()) % 2 == 0:  
#        return (srtd[mid-1] + srtd[mid]) / 2.0
#    else:
#        return srtd[mid]
    
    fig = pylab.figure()
    ax = pylab.gca()
#    a = sorted(S_phase_timings.values(), key=len(vals) in S_phase_timings.items())
    bp_S = ax.boxplot(S_phase_timings.values(), 0,'gD', patch_artist = True)
    for box in bp_S['boxes']:
        box.set(color="#000100", linewidth=2, facecolor="#FFFA03")
    for whisker in bp_S['whiskers']:
        whisker.set(color="#000100", linewidth=2)
    for cap in bp_S['caps']:
        cap.set(color="#000100", linewidth=2)
    for median in bp_S['medians']:
        median.set(color="#0002FB", linewidth=2)
    for flier in bp_S['fliers']:
        flier.set(marker='o',color="#0001FF", alpha=0.75)
    ax.set_xticklabels(["(%s, %s)" % (len(vals), t1) for (t1), vals in S_phase_timings.items()], rotation=90)
    ax.set_title("S duration")
    ax.set_xlabel('Positions')
    ax.set_ylabel('Duration (h)')
    ax.set_ylim (ymax = 6, ymin = 0)
    pylab.tight_layout()
    fig.savefig(pm.output("S_duration.pdf"))

    
    fig = pylab.figure()
    ax = pylab.gca()
    bp_G2 = ax.boxplot(G2_phase_timings.values(), 0,'gD', patch_artist = True)
    for box in bp_G2['boxes']:
        box.set(color="#000100", linewidth=2, facecolor="#71FA03")
    for whisker in bp_G2['whiskers']:
        whisker.set(color="#000100", linewidth=2)
    for cap in bp_G2['caps']:
        cap.set(color="#000100", linewidth=2)
    for median in bp_G2['medians']:
        median.set(color="#FF0104", linewidth=2)
    for flier in bp_G2['fliers']:
        flier.set(marker='o',color="#0001FF", alpha=0.75)
    ax.set_xticklabels(["(%s, %s)" % (len(vals), t1) for (t1), vals in G2_phase_timings.items()], rotation=90)
    ax.set_title("G2 duration")
    ax.set_xlabel('Positions')
    ax.set_ylabel('Duration (h)')
    ax.set_ylim (ymax = 6, ymin = 0)
    pylab.tight_layout()
    fig.savefig(pm.output("G2_duration.pdf"))

    fig = pylab.figure()
    ax = pylab.gca()
    bp_G1 = ax.boxplot(mito_unknown_timings.values(), 0,'gD', patch_artist = True)
    for box in bp_G1['boxes']:
        box.set(color="#000100", linewidth=2, facecolor="#009704")
    for whisker in bp_G1['whiskers']:
        whisker.set(color="#000100", linewidth=2)
    for cap in bp_G1['caps']:
        cap.set(color="#000100", linewidth=2)
    for median in bp_G1['medians']:
        median.set(color="#FF0104", linewidth=2)
    for flier in bp_G1['fliers']:
        flier.set(marker='o',color="#0001FF", alpha=0.75) 
    ax.set_xticklabels(["(%s, %s)" % (len(vals), t1) for (t1), vals in mito_unknown_timings.items()], rotation=90)
    ax.set_title("G1 duration")
    ax.set_xlabel('Positions')
    ax.set_ylabel('Duration (h)')
    ax.set_ylim (ymax = 40, ymin = 0)
    pylab.tight_layout()
    fig.savefig(pm.output("G1_duration.pdf"))
    
    
    fig = pylab.figure()
    ax = pylab.gca()
    bp_M = ax.boxplot(mito_timings.values(), 0,'gD', patch_artist = True)
    for box in bp_M['boxes']:
        box.set(color="#000100", linewidth=2, facecolor="#FF00FF")
    for whisker in bp_M['whiskers']:
        whisker.set(color="#000100", linewidth=2)
    for cap in bp_M['caps']:
        cap.set(color="#000100", linewidth=2)
    for median in bp_M['medians']:
        median.set(color="#FFF111", linewidth=2)
    for flier in bp_M['fliers']:
        flier.set(marker='o',color="#0001FF", alpha=0.75)
    ax.set_xticklabels(["(%s, %s)" % (len(vals), t1) for (t1), vals in mito_timings.items()], rotation=90)
    ax.set_title("Mitotic duration")
    ax.set_xlabel('Positions')
    ax.set_ylabel('Duration (min)')
    ax.set_ylim (ymax = 100, ymin = 0)
    pylab.tight_layout()
    fig.savefig(pm.output("mito_duration.pdf"))
    
#    with open(pm.output('mito_duration.txt'), 'w') as file_txt:
#        wr = csv.writer(file_txt, delimiter='\t')
#        wr.writerow(mito_timings.keys())
#        wr.writerows(izip_longest(*mito_timings.values(), fillvalue=""))

   
    
    

   
    
    fig = pylab.figure()
    ax = pylab.gca()
#    sort_try = sorted(mito2mito_timings.values(), key=(vals for vals in mito2mito_timings.items()))
    bp_Cycle = ax.boxplot(mito2mito_timings.values(), 0,'gD', patch_artist = True)
    for box in bp_Cycle['boxes']:
        box.set(color="#000100", linewidth=2, facecolor="#FFFBF8")
    for whisker in bp_Cycle['whiskers']:
        whisker.set(color="#000100", linewidth=2)
    for cap in bp_Cycle['caps']:
        cap.set(color="#000100", linewidth=2)
    for median in bp_Cycle['medians']:
        median.set(color="#000100", linewidth=2)
    for flier in bp_Cycle['fliers']:
        flier.set(marker='o',color="#0001FF", alpha=0.75)
    ax.set_xticklabels(["%s_%02d (%d)\n(%s, %s)" % (w, p, len(vals), t1, t2) for (w, p, t1, t2), vals in mito2mito_timings.items()], rotation=90)
    ax.set_title("Cell-cycle duration (Meta to Meta)")
    ax.set_xlabel('Positions')
    ax.set_ylabel('Duration (h)')
    pylab.tight_layout()
    fig.savefig(pm.output("mito_2_mito_duration.pdf"))
    
#    fig = pylab.figure()
#    ax = pylab.gca()
#    ax.boxplot(mito_unknown_timings.values())
#    ax.set_xticklabels(["%s_%02d (%d)\n(%s, %s)" % (w, p, len(vals), t1, t2) for (w, p, t1, t2), vals in mito_unknown_timings.items()], rotation=90)
#    ax.set_title("G1 phase duration")
#    ax.set_xlabel('Positions')
#    ax.set_ylabel('Duration (h)')
#    pylab.tight_layout()
#    fig.savefig(pm.output("mito_unknown_duration.pdf"))
    
    pylab.show()
    
#    ax.scatter(mito_timings.values())
#    ax.set_xticklabels(["%s_%02d (%d)\n(%s, %s)" % (w, p, len(vals), t1, t2) for (w, p, t1, t2), vals in mito_timings.items()])
#    ax.set_title("Mitotic duration (Meta and Ana)")
#    ax.set_xlabel('Position')
#    ax.set_ylabel('Duration (h)')
    


    
    print 'Fini'
    