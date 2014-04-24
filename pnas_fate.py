import numpy
import pylab
import pandas
import matplotlib
import cellh5
import cellh5_analysis

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
    
    def combine_classifiers(self, output_name):
        all_combined_classes = []
        for _, (plate_name, w, p, track_ids, track_labels) in self.mapping[['Plate', 
                                                                            'Well', 
                                                                            'Site', 
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
    print 'Start'
    plate_name = '140414'
    ch5_main_file = "/Volumes/groups/gerlich/members/Igor Gak/FINAL/analyzed/hdf5/_all_positions.ch5"
    pos_mapping_file = "/Volumes/Seagate Backup Plus Drive/Screening data/FINAL/FINAL.txt"
    hmm_xml_constraint_file = "graph_igor_pnas_h2b_combi.xml"
    
    pm = IgorsClass('igor_test', 
                        {plate_name: pos_mapping_file}, 
                        {plate_name: ch5_main_file}, 
                        sites=(1,),
#                         rows=("C", ), 
#                         cols=(2,),
                    )
    pm.read_events(5, 99999)
    pm.track_full_events()
    pm.combine_classifiers("Event labels combined")
    pm.setup_hmm(5, 4, hmm_xml_constraint_file)
    pm.predict_hmm("Event labels combined")
    
    pm.plot_track_order_map(["Event track labels", "Event labels combined", 'Event HMM track labels'], 
                            [IgorsClass.cmap3, IgorsClass.cmap4, IgorsClass.cmap5])
    
    print 'Fini'
    
    
    
    


