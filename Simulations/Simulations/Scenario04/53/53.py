#!/usr/bin/env python
from numpy import *
data_1 = [
array([30.946033390341114, 31.146125942680833, 31.870932309113833, 32.20905215901647, 32.60075915069947, 33.709584926179225, 30.89800227677789, 28.26769980068224, 32.91078752835131, 33.15584168423043, 32.04695563358184, 30.923252665347505, 31.568884312226626, 34.281301124215716, 31.755031277445443, 31.77340218001824, 33.41074570119948, 31.801616654777323, 31.315849860402313, 34.293437977493795, 31.076718957051614, 29.70435764654914, 32.27552231099221, 32.23228877687726, 33.66421456055807, 29.261346506647207, 33.35360541230461, 33.06268003282788, 28.9896787106629, 32.653879414546914, 28.216733740511533, 28.45262498466057, 33.00151077955063, 31.122936934685132, 31.720413136470075, 31.528731577827728, 33.2808186139784, 29.562121564998375, 31.037522761395966, 32.82596114338915, 31.8303500987488, 33.62149836543714, 31.11579326934082, 25.984875753134776, 24.8244729005064, 24.73879595671341, 26.88879530277641, 25.356648615109897, 27.01778463821762, 24.854535379324222, 24.472506994090086, 25.811420310648426, 28.0223056864666, 23.593133972561844, 26.19986256813676, 25.68194799150365, 23.885726348125424, 24.46730708314181, 22.17627326220955, 22.663525801934007, 22.640603636427116]),
array([33.338772311162536, 33.404275617172985, 32.28282097202031, 34.1849244030094, 33.56203092978748, 33.227191082399656, 33.51327024873906, 33.38644509212199, 32.377143984134754, 34.8113027389262]),
array([31.839753783849076, 30.582126517538065, 30.833645958869877, 32.3095648294035, 31.350442238000507, 31.74399327806097, 29.415466437040994, 30.614311763215355, 31.576371894790075, 31.941798468153, 32.149509270292384, 28.964487783816406, 30.987653078546227, 31.78426603957817, 30.31074100422485, 30.63889142958109, 30.80638329657037, 30.560151464400434, 31.53600135572605, 29.44007267300305, 30.68822549914273]),
array([28.41132610213543, 29.735379550661634, 28.91977634105861, 28.525536663340496, 29.35524715871085, 28.296994904840346, 29.41771807998585, 30.570934242374346, 28.19653968832109, 28.41837549081609, 29.306455069836, 28.400769483629315, 29.38337813765319, 28.745241771077747, 29.25811243436631, 30.437225533839904, 29.932242146508525, 28.694128518270354, 28.686884786866223, 28.7105791006213, 28.955182879144495, 28.907236292127195, 30.040405262791353, 30.466388229670805, 28.401299707285975, 30.16045001544684, 29.309946910729128, 29.753280908030096, 29.112761168146886, 28.64041289677595, 28.512656604163407, 28.162710794285964, 30.15426263375966, 29.825739340561594, 28.777002462870716, 30.228052481232847, 29.91326799891592]),
array([29.56715730301766, 29.245757094621556]),
array([28.6422938878021, 28.263042839231904, 28.601624549243144]),
array([28.9449738240909, 28.670541917975914, 29.516995898870043, 29.108907667641553, 26.92167349541556, 24.380353668488176, 25.126486594762444, 26.82498905206475, 25.57447469338634, 25.419731065500024, 25.892575021085953, 27.223264255473524, 21.78411698063405, 17.569550759700032, 20.250649407348497, 15.267539447874501, 15.81690939229309, 12.261824948758626, 11.720132919413878, 12.14398228318429, 11.865022872605852, 11.386152561995008, 8.178226402349171, 10.531535225368561, 9.508005605034056, 1.940300759648857]),
array([26.970391623748032]),
array([28.18779793564632, 28.29543770058036, 28.160327675867975, 27.778834729132946, 24.750504787100347, 27.685248827708012, 26.31299107906821, 26.83900534855273, 26.844738829959205, 23.22132221986063, 24.072879367686188, 27.69074886525829, 27.025500032587754, 26.073417386839928, 23.08849775114185, 26.878279006627373, 25.726483872724916, 26.631879526726234, 21.009196401432575, 20.477469791735192, 21.08876177745705, 21.19762267716181, 21.60227922032435, 16.849503054207293, 18.877744507342108, 18.168041623458148, 19.636401674350502, 18.96500634977831, 20.251388993386712, 20.32132087223746, 20.40087916258253, 17.9192057440183, 14.095537642213966, 15.25241811591091, 15.576119208417724, 14.271355103715583, 15.247016392163411, 14.789798226894579, 15.218843792400646, 14.55476543140935, 12.20304859818523, 12.733022503056848, 12.325685276721535, 13.243349234535637, 12.377434491229037, 9.326626592896677, 11.46290608642481, 8.857392911260249, 9.043546424693176, 9.537128783853372, 9.332183434223468, 10.96991324736193, 9.241919749099248, 7.791105331836395, 9.687180397695785, 10.458242332896397, 10.830328930396364, 10.124355715004045, 11.461404189975081, 10.438013848953418, 9.549302050337637, 9.234221739746143, 9.091713091029447, 8.470239080288946, 9.739320082309234, 11.277331267507018, 11.011060602889872, 7.8749836079244515, 7.726441111144874, 8.229986970076117, 8.119818332564606, 10.015109936809715, 9.33589697270543, 9.143310910520793, 11.611661497788752, 6.352519984531686, 6.521056657778113, 6.9184777092789105, 5.627865244397318, 4.0955434900344425, 3.879362125527005, 5.124274061836515, 4.3454160327973845, 4.261271496732739]),
array([26.659612134476717, 23.073432433930684, 25.24364779180022, 22.941570772366713, 16.79608730753574]),
array([25.51691664449089, 24.740841908683052, 25.346777247624935, 26.20010700588456, 25.086370437247503, 23.519738806324817, 26.65742743520682, 23.184075374290444, 23.143408878223955, 25.29322849220295, 26.31607955351841, 22.84290468826912]),
array([22.079786339721515, 17.32144828480003, 17.39756317942011, 15.859787477622087, 15.823776213681889]),
array([24.066500360056533, 23.374229400059267, 24.28908738241156, 23.19683328916708, 24.05359960513713, 22.888906304315295]),
array([24.13650669495088, 23.805794874080135, 23.9010434868401, 23.211514424654055, 23.95993836158972]),
array([23.292455822605096, 22.719833815927714, 22.9067305333899]),
array([23.040555914173773, 23.118835419354983, 22.345765212662272, 22.93719794476713, 20.95329627379407, 21.19503665005297, 20.15610459441917, 17.260585867523957, 18.402905252645294, 18.740567460601735, 19.833479988668607, 18.77072680570472, 15.83995030934012, 13.940164842356221, 12.414286710045706, 11.657275589985984, 9.386894174616678, 7.788807195380424, 8.473796284491177, 9.386414651439589, 9.065415931409463, 11.084857959142491, 10.516466139461064, 9.626327459043676, 10.803094090009184, 10.44111038425364, 9.072557852876846, 10.294939847233536, 8.109933168062598, 11.624431180300391, 9.87516224632722, 11.102960788770124, 9.257271094289615, 7.911657578272356, 11.090819913372657, 11.181769583169007, 6.544153007460365]),
array([21.912852873096323, 22.855354593145528, 21.377150946990476, 21.03379890126161, 22.936792036578638, 21.82772089556883, 21.732251602925466, 20.672054048812242, 22.833428004408884, 22.658037213631363, 18.98241761153119, 18.024373662990364, 18.277709471509727, 18.038104047863058, 19.62077960168716]),
array([22.362963549281044, 22.074577573194034, 22.17519796848364, 19.492628031862417]),
array([22.43412882540237]),
array([20.64557291899419, 21.211995308322354, 20.970548103918876, 21.253539315288567]),
array([21.14052278845559, 18.279190994898357, 18.33198008866313, 14.138029130786329, 14.336791412820597, 12.058161660029754, 9.123813890223293, 6.552340585831944, 6.101728886783754]),
array([20.00937446979228, 17.526020586689334, 17.33006109627036, 17.60619779308399, 17.981074759690756, 20.369524230468084, 19.411061577321423, 17.695299631281546, 15.51101740520357, 15.893190770955083, 15.340392959159566, 15.806028362300587]),
array([18.63351800698856, 18.415715343404592]),
array([17.212333491808515, 18.27262379364684, 18.51495737990962, 16.224761659733126, 16.487259452640863, 18.11168641390359, 17.20756859166905, 14.313407584106939, 15.361073986496299, 14.522124989580899, 15.157689839525235, 15.17710649177394, 14.246598342064452, 13.068977093730751, 12.914445362447436, 11.634536442603984, 12.68341560553377, 13.675412480653204, 12.506822633361137, 11.032417047434532, 11.332441816093423, 11.177419798991313, 11.237168309330212]),
array([18.77560440873397, 17.380172969480437, 17.648085997951156, 14.499863370738467, 14.296334144046508, 14.330261317059385, 13.384915922258637, 10.119059244680185, 10.80324052820735, 7.543854302371617, 8.377712733582996, 10.668611078146752, 9.364727904522411, 8.356204395748694, 10.51861164016212]),
array([16.867529786620906, 16.28829853311729, 17.344103752702463, 15.973677911805323, 16.476682539416814, 16.90748843611978, 15.48209988142448, 14.793976481411082, 15.61604723551445, 14.295631557789575, 11.63692110131053, 13.57286006408859, 13.135375094715538, 13.402272243240564, 9.379500422140582, 8.351361725362194, 10.091495150585988, 9.579438331928435, 10.534451918378377, 9.854969788040313, 7.9885422316302, 10.264867211129214, 7.390389431060943, 9.314179405348177, 10.668569489415454, 11.170066587340168, 8.001134159671835, 7.526439831207986, 9.260117206191033, 6.38491574999687, 6.590773156762458, 6.264267964540197, 5.031107338177603]),
array([17.51428466633977, 16.5118344170367, 18.265716250795847, 18.379821268372304, 16.760790034966856, 17.67261615378979, 14.883212810581213, 15.603523846785116, 15.467656869011057, 14.088778583517122, 13.657397463859246, 12.594121640182232, 10.512584066383285, 8.438576743956895, 8.319839039162561, 10.384695830609955, 8.441491324729299, 8.693616401025556, 8.867727527527475, 8.82423834455639, 10.983563378087728, 9.939300677848884, 10.678461217014009, 8.658771906078861, 11.52427130930791, 9.228012376126099]),
array([7.584554064338713, 9.56962905156533, 8.879915946557256, 7.079970202336626]),
array([16.014643670318918, 16.701153412354067, 16.716807335039253, 17.001230056742582, 14.409332935125988, 14.249334801070065, 15.587140222664377, 15.466735081396141, 14.62883799410191, 15.138173450939187, 14.362923278748415, 11.681513140888377, 11.794255710190498, 12.786696595136672, 11.918620114714454, 13.509809573464985, 11.928859045806743, 12.829294587861359, 10.957075692891976, 11.390924961812852, 11.344180352591241, 10.773853561319738, 11.154969144517644]),
array([16.177389031990412, 14.946749936311038, 15.628341940499785, 14.813210797807198, 14.004780066523715, 15.259764329673454, 14.355240082285398, 14.327865859812393, 15.261469782554371, 12.649633361076205, 13.651664579171161, 10.285678070370302, 9.86809067253266, 9.841189949824916, 10.718950667469466, 8.556356686850833, 10.444811340540229, 7.550920748443571, 7.85577934664075, 9.718906724441592, 10.235456726072453, 7.535632028670432, 7.808619509012205, 10.644054926114885, 8.582757427533549, 10.322141740445662, 9.150957956539544, 8.804081775447276, 10.093265383592287, 10.884491311558444, 9.00777265791403, 9.015791207278038, 11.267325914968996, 10.67400878306257, 6.814326327452397, 4.722404092730777, 4.233725327445892, 4.82505801122828, 2.590088938719422, 3.354821172723997, 3.535571717952652, 3.515906220024399, 1.987925246958746, 0.7903824211167096, 0.0]),
array([15.008303458754805, 15.015219229095672, 14.694175926220817]),
array([16.276102428014745]),
array([15.16389858440367, 13.922777136041526]),
array([14.661013893652827, 14.717714155259877, 14.603834371687658]),
array([14.15198081601153]),
array([14.869396243399105, 15.174767315986552, 14.518078908744343, 15.190262644616489]),
array([15.060809194679619]),
array([14.819049050642846, 14.958278742397228, 14.639269214016341, 14.456490230616303, 14.880596904199793, 14.659422081651709, 14.930130276318465, 14.874387954568432, 14.713315945004936, 14.430366105857287, 14.912439286668777]),
array([13.830443976884522, 14.36575900136251, 14.54550140909491, 14.007867231024504, 12.06309109119586, 13.20064997762984, 11.75024480759419, 13.50124742381638, 13.504741414315525, 13.099753148705094, 13.7545771516793, 12.073827418684745, 13.528832326962306, 13.506912848618313, 12.319086889302882, 11.084520628201068, 10.327369709939118, 10.703571713962006, 10.285094631174326, 10.130716826263615, 11.599402631747616, 10.26310969715121, 9.775917757895797, 9.335129738536207, 10.950154047785938, 11.202695814199, 10.007198866746283, 8.934805501828755, 9.226927141148796, 10.14715099962864, 9.323826683787688, 8.908668961767809, 11.279711655365816, 9.434541307432681, 10.733940292404629, 9.3888731870913, 11.406992360235082, 10.996762892280614, 11.166825049303911, 10.462564357340344, 9.5345988087794, 8.687522321141461]),
array([14.180720660104814, 14.121222269345832, 13.281527358150964, 13.763934611273706, 12.693856735013597, 12.254286085586337, 12.631556975147925, 13.078953493814257, 13.443112162430271, 9.345313884135308, 7.7145594002098345, 8.912881106819375, 8.698715107926668, 8.086447665296024, 10.826752731476704, 9.357382291527475, 11.366352557788987, 8.791987648489084, 9.86488717488138, 8.493964871962142, 10.333142607859315, 7.127173702893209, 7.1647563064882345, 6.73056342965369]),
array([13.86990569922662, 13.99682934272443, 13.72624697347226, 11.761724222860117, 12.441531136273708, 13.741338572927022, 12.314791577575834, 12.140570146637936, 12.42564465200462, 8.503750398617967, 9.458317481555804, 11.053410439383656, 8.596601383446307, 10.055035043383445, 10.864298280239929, 10.250971050659208, 9.484675963515363, 9.971105281125274, 8.909441900941186, 10.299397669839164, 8.599370701428665]),
array([12.546809769659836, 11.394183335352588, 11.585207581655007, 11.544360380111456, 8.674097764424186, 11.548701039827893]),
array([8.457665599318473, 11.467617089963532, 10.262189761189685, 10.72056180951088, 8.049513519163751, 10.605663136097437, 8.01567344211766, 10.172624680269555, 11.227138736907404, 7.956576971511316, 8.708305113403323, 9.31793401804648, 8.324082519566042, 9.700031023013402]),
array([13.033549992005081, 12.537230806619775, 13.071173581700425, 11.231276691780005, 9.73077474009328, 11.095052140636758, 10.446823744691066, 11.262860495062395, 10.775702147489959, 9.223913838527983, 11.595632029598017, 10.942564745969962, 10.002405302691782, 11.408511951499447, 10.52139385499956, 11.605418056072383, 10.509293565709363, 9.694772086284743, 9.76546989077191, 11.44319980594831, 9.556135342005785, 9.476105343441592]),
array([12.574542801402933, 12.450678429914813, 12.043642040495609, 12.716242729667105, 11.653439547531057, 12.785889495168467, 7.355796076058472, 8.615014786380375, 8.487742889042636, 11.020345061118318, 8.614743298532636, 11.05451511428452, 8.267721919170324, 7.491244147134781, 10.72519010784113, 11.015944123299905, 8.710237534627675, 7.2834453480424575, 10.70000774073171, 9.83050346430102, 7.727439206054868, 7.819554337262542, 10.359926769502763, 9.710451989233974, 10.158399177778659, 7.890823675936183, 10.277082883948033, 9.729450722899976, 7.508782382648339, 11.455952180338809, 10.538444304786255, 7.833307151366919, 7.811258594314153, 11.515021556069815, 10.506144174415851, 8.989163708016754, 11.272419262734546, 11.468044742196676, 8.965220403890081, 8.833543359605512, 6.307068012342689, 4.154691398318651, 4.514832659429427, 4.66231667899132, 3.4019821764605624, 3.54412988022894, 3.569216743496307, 1.2856776218952761, 0.0]),
array([12.487848968391518, 12.350470917612405, 12.111095563032462, 12.550633702025795, 12.221444218252563, 12.048344321273937, 12.306647687835868, 11.674249400104213, 11.35151828512052, 11.350852911900827, 11.278031595690186, 11.486745847133992, 11.111837403161868, 11.273503748720767, 11.254939434448637]),
array([12.094450701351557, 11.656943628631351, 12.315430178315781, 10.853887163287382, 11.587966293580688, 10.992299753630643, 10.610322160929613, 11.078274115064525, 10.80283419516985, 11.61128160347404, 10.418151515377645, 10.6107631305612, 11.080027457508304, 10.31299142018948, 10.712899497379066, 10.285231090125075, 11.321646435291198, 10.397943316385897]),
array([12.243416767617111, 11.655019435148077, 12.17226225920448, 12.027145324416363, 11.751238965324365, 11.888356093845536, 8.990744065729071, 10.352068042650792, 7.9623401473748165, 8.261546571866443, 8.72928470303139, 7.931098878035221, 9.618711119221985, 7.561084720296223, 8.264800959785338, 8.40062904026388, 8.89173571662306, 10.337103429685177, 11.156277738700712, 8.828539732015374, 8.451730363376496, 9.376426601713309, 10.472120149649994, 8.865276687782503, 8.817716137710828, 7.317911898079837, 11.620534424602976, 10.464002765465198, 10.95742818040673, 9.61517781320023, 10.865885398576008, 7.405090222019351, 10.736109380861015, 9.04782425149147, 8.416180433922573, 9.502567335111461, 9.45963904774138, 7.865352083696659, 8.983358413575312, 10.841755658405708, 7.69614705320182, 10.183426487768218, 7.325228369696618, 8.097335768078644, 6.506396741737077, 6.783849858320282, 6.901087560459459, 6.6060642341604945, 6.772369730515332]),
array([11.767922554779716, 11.75317061836715, 11.59629846100778, 9.238662414365148, 9.32935568233019, 10.969056504450409, 7.726785734274912, 10.22673382476748, 10.335403003052964, 7.585804230039748, 11.389666901236403, 9.95766837131481, 10.464348830434249, 7.802211819485466, 10.65121247566259, 11.02117516526747, 11.451312418382424, 9.314740053345767, 8.85422037128906, 8.157583497554896, 11.586124766115056, 10.43164833050108, 9.001193255984202, 11.52623573688993, 10.345926099586743, 7.685860918130468, 8.235432726539909, 7.859888322581703, 7.941797422904212, 9.567147641150264, 9.43115778134562, 11.173299792593257, 7.40249143898921, 9.109574341662944, 7.639187903480691, 11.062498428467174, 8.027298259164109, 9.119117757531221, 9.211075837236699, 8.514265593559728, 7.520236813259416, 6.681094128826008]),
array([8.57609700910492, 11.5446227593403, 7.868886670058785, 9.838421414723923, 7.830514564640004, 8.049166564916305, 10.347551891939082, 10.596652948781013, 8.404164873687774, 9.19027129495651, 10.10868025877772, 10.70548955352287, 11.539403453978018, 10.694043376903078, 11.409636363848193, 9.389468127341495, 8.405982151256786, 8.09455724551396, 6.2423427474474655]),
array([11.43200804559219]),
array([10.474929519393136, 9.504613603680076, 9.777362457182486, 10.040066971132138, 10.158772362272634, 9.699032003129691, 11.160584643948978, 9.971222543690807, 10.937513586803167, 10.36876947236858, 10.236824917729994, 10.863679063070729, 9.984579636396058]),
array([10.690709562789786, 8.136213314763078, 7.855235164042703, 9.883510450124435, 9.824039881743701, 9.15932291858364, 7.967468927315656, 10.723271412338393, 3.9756731641474925, 2.004256824587355, 0.8855878106604402, 1.4530562846475692, 0.0]),
array([7.070507796117607, 1.5333956426206383, 0.0]),
array([9.440870705227471, 8.021216559957491, 7.387523390661119, 8.312654629430337, 8.966695480604368, 7.588467726105989, 7.95229753225904, 8.663746378973437, 3.3427024898379156, 0.0]),
array([8.899694764111715, 9.957395317560946, 8.570910639109997]),
array([8.698771508131095, 9.185753331045133, 8.78156257151317, 8.32518463646607, 9.917882960345827]),
array([7.284905295390471, 8.419401897501912, 9.822252943956931, 8.923251596658625, 7.525704123105829, 9.670401265548728, 7.593859566333785, 6.728922089861045]),
array([3.8307673023831414, 3.556832224703862, 0.0]),
array([8.358737036927177, 8.93644232421279, 8.23917325499718, 8.734193051758766, 8.261322237249571, 7.482013942783981, 8.806054881855228, 8.707898648355593, 3.9416081503582467, 3.808820776981463, 5.079359065752559, 4.323166352113966, 4.566124676461324, 2.5512505286612956, 0.11247392947055038, 0.0]),
array([7.335577359703326]),
array([8.763673094075944, 7.505375858644433, 8.321112510897287, 9.167283320476187, 7.795973618451373, 9.248378475454128, 8.288723785643077, 8.699751072530653, 7.170717953179281, 5.990342319659881, 3.614208899167766, 4.087463727107542, 1.9469114645816061, 0.12276795187183571, 0.0]),
array([7.516456242662856, 8.168593396930314, 9.374940585105461, 9.31084813054409, 9.1623380061887, 8.986939716579897, 8.078116880894902, 9.315806570504819, 7.295993013991439, 7.251378840119338, 7.868475907799218, 8.248522158702066, 7.868129718228049, 9.18339626512248, 7.896478008209783]),
array([9.078337643355637, 7.737679509071069, 7.523133539897343, 8.711083541997027, 7.31680228867044, 7.471700871792206, 7.039027662435983]),
array([8.226035743699839, 7.419517372041067, 8.535625240127864, 9.041687427020214, 7.741760040479442, 9.063244346965508, 5.853799864972269, 6.709678860943458, 7.091502710762704, 5.101027150329607]),
array([8.8144477752073, 7.6224567796652405, 8.366650517082421, 9.081928604695998, 8.064366921125853, 8.274094862943699, 7.8151400655787615, 8.893253547044731, 3.9085486577935686]),
array([7.312442534456844, 8.525190978885972, 7.2997864330332645, 7.382896647042999, 7.270178708310807]),
array([7.9893517223210155, 8.633959962535833, 8.706615350957678, 6.982974451318907, 3.640830786818567, 4.087867483049344, 5.024158987728382, 3.1345917347029895, 3.420495246073428, 1.853793651538431, 2.2050837059296926, 1.1496448399308452, 1.7678095104447182, 0.2125614269528422, 0.06595083059337584, 0.0]),
array([8.994079923493063]),
array([7.945052187193668, 8.21753249191463, 9.166604478038867, 8.62582231941561, 9.216285794699544]),
array([7.266534083101106, 7.88630325395446, 8.944116732441945, 8.165679933706695, 7.749337006875088, 8.89115998621592, 8.137955923072635, 8.958946013719796, 7.31903419901891, 8.779013105277198, 8.658935968103174, 8.030575570243588, 9.13028962250355, 8.865586534505958, 8.545356308971705, 8.22430792244238, 5.389703170815626, 6.035188351421249, 6.52001210146912, 6.978143972591869, 6.655161892081305, 7.220494462092991]),
array([8.11442970221487, 0.0]),
array([8.81486656346295, 8.854054226231481, 7.65736965407885, 7.404639488563591, 7.361895793705011, 7.572821794639624, 8.4059871056565, 8.58055627077024, 7.841127237641241, 7.426060306682608, 7.412943862301324, 8.908875244637928, 8.072160321529427, 7.93428488738944, 7.456182645822991, 7.518340801283831, 7.633494895761616, 8.212851768433998, 8.231332796074252, 6.6622818111120115, 4.8982752243328935]),
array([7.610278635496988, 7.525997433196963, 6.765260166437396]),
array([8.398400284404447, 7.798856846550079, 8.030290673082643, 8.125075979260567, 7.407301902081715, 8.159527105831547, 7.7913733890802, 7.3997171153328125, 7.252381800778211, 8.114850896464786, 7.675201297525684, 8.546555399523806, 7.725699965551367, 8.472086486501778, 7.774467642049934, 7.427696053650351, 5.880595382731084, 5.706025387623136, 6.195101411849027, 6.536081900213781, 6.750371373045126, 6.369799737094172, 5.498701927436379, 5.832099325894912, 6.874697126685691, 5.237494113621519]),
array([7.9002489197369705, 8.07824927998914, 7.883112744668841, 7.948375714310572]),
array([7.282618835511538, 7.744888769848504, 7.114302416352166, 3.473740007597298, 2.275980756783621]),
array([7.815215757652291, 7.754324257357873, 7.408541121620157]),
array([7.2744481449492095, 7.5006138620954665, 6.202557957039837, 5.953001734605735, 5.495338171859965, 5.806225232203189, 6.346070881128801, 2.877600179400114, 2.990448814001768, 2.2666294139799885, 1.8204780366979916]),
array([7.360452066928549, 4.214878746474659, 3.631104361428095, 1.8271436198780333, 0.0]),
array([5.093735816041668, 5.117405769419629]),
array([5.576560823331915, 5.603682110200926, 5.52656410947235, 6.643379042404692, 5.07537449452, 3.8758057488675512, 4.013835752216289, 4.543260293850075, 4.668135384766052, 3.651407564103558, 2.847456875022778, 2.752648913418673, 2.4777005108901036, 1.6765173672876497, 0.16124285383924752, 0.0]),
array([2.9198272575066886]),
array([3.713854364304052, 3.579153997320225, 0.0]),
array([6.07406294878469, 5.658962541808289]),
array([6.627793990466143, 5.147851774604312, 3.9764899825987468, 4.828897594060024, 3.1422812153124755, 2.71112894863453, 3.2121227609649847, 2.623372405352167, 1.7572378364010521, 0.4171156514141, 0.1041370234383773, 0.0]),
array([6.443451467254533, 6.655742291556594, 5.0663688067547685, 4.2394197228825155, 2.057633191218326, 1.7489168222460036, 0.3731843222312701, 0.05933899621075042, 0.0]),
array([5.748267041932992]),
array([6.70347495608714]),
array([6.118736639034006, 5.542896680706586, 6.656043176525692, 5.767761496321534]),
array([1.2317479204565118, 0.4891039675972665, 0.0]),
array([3.652407375898826, 0.0]),
array([6.581458485473356, 6.126155079203197, 4.769047437142113, 4.769653064739776, 3.92382733974128, 2.7889747865403853, 0.27528502913016584, 0.0]),
array([3.637101139162989, 3.214848213539455, 1.7516318263608028, 0.0]),
array([5.661110405003491, 4.261087666835561, 5.309024963227534, 2.318329755443573, 1.088708246340329, 1.214572915865563, 1.7044234923237074, 0.0]),
array([6.017449575009865, 6.081447622691149, 3.6673489841384965, 0.0]),
array([4.524136071227911, 4.053363272860323, 5.162365909935298]),
array([6.287806665993492]),
array([5.34403302852973, 3.812161360761568, 4.187071014277327, 1.2068852736081506, 0.0]),
array([4.520675947432439, 1.8134379390817768, 0.0]),
array([5.0474004105822035, 0.0]),
array([5.698552280467227, 4.51532582085974, 2.0355260442339658, 2.5135713025561057, 0.0]),
array([5.356190568552403, 3.665557905817394, 4.3293413193316965, 5.276027255304342, 3.2155512430806192]),
array([3.3833368786119706, 3.3571032117139654, 1.4368645279010488, 1.0136651267825911, 1.6829019729357477, 0.05862948673585876, 0.0]),
array([4.539281246372161, 1.2906093498350009, 0.01713250408395872, 0.0]),
array([3.855193712870717, 4.615693558372659, 2.886614746970961, 1.856953639997548, 0.9764180313884692, 0.0]),
array([3.2947828283495637, 3.392141489554763, 0.0]),
array([4.065847646366367, 0.0]),
array([4.606511080650324, 4.214352738986348, 4.245958344017825]),
array([0.6197015088125959, 0.0]),
array([4.525409048108281, 3.9671898613961147, 3.6998066511642946, 3.7504262564130952, 2.755746232578668, 2.64547954207103, 1.946994974623272, 1.7127473284861894, 1.483438063493212, 0.9361729489563496, 1.31822675535483, 0.7012994773966619, 0.5909662944367198, 0.2738732635413209, 0.34172663233131684, 0.1003775306689073, 0.0]),
array([1.2078649695897594, 0.0]),
array([4.060350049521153]),
array([3.4088422764748474]),
array([4.401829196111464, 4.3715337446274285, 3.8822989339047016, 3.964098072331275, 3.6888674802090025, 3.156379628879004, 3.3445275517886186, 3.173708285570786, 3.471793166945269]),
array([4.683425591758137, 4.54599662928164, 3.6027826338423576, 3.913817468980798, 4.052525736152268, 4.250246539989822, 3.966651279195774, 3.50762745018832, 3.5816433705988424, 1.047353018776473, 1.0425441686352426, 0.9297105531147997, 1.7987732746931198, 0.7893867945455142, 1.5682672063527576, 0.5564488804808072, 0.33178433376287203, 0.0]),
array([4.6479994404872205, 3.6104646900189277, 3.735576880797683, 3.8418769420513685, 4.07635577380401, 2.7593148913587413, 3.534912004265996, 2.0135612942919385, 1.3214091547874554, 0.39947323680681285, 0.0]),
array([3.707139730519475, 4.335381042326007, 4.45514609596427, 4.213992171735565, 3.7132136791259365, 2.94848339837002, 2.8544547826750204, 3.4439986309469983, 2.7049034677976653, 3.087430625577875, 2.671105860979191, 2.4162080541802737, 1.4047353494437627, 1.2185686125728201, 1.6950076012893196]),
array([3.7946459945276105, 3.93139138682357, 4.35069992139584, 3.8674301290040902, 2.821215713796457, 2.6649256645684973, 0.9590980540078158, 0.0]),
array([3.7936557294830577, 3.6582120869118544, 0.0]),
array([4.119442178785926, 3.415358534493774, 3.3925879862869106, 2.0456585231346898, 1.3483468886338532]),
array([3.88794157678197, 3.496883252366122, 2.5103931870260268, 0.09517946127748748, 0.0]),
array([3.6201619966347756, 2.658995081259682, 3.473748122937098, 3.2730817813625404, 2.5474086161073712, 1.84003361830847, 1.0201061189411673, 1.4873632496176645, 0.6601595337589355, 0.0]),
array([2.60692276110317, 2.499654982570752]),
array([3.0035866323062996, 2.545652443479503, 2.4841868164525733, 0.0]),
array([2.791820954790767, 1.037785740754582, 0.0]),
array([1.352273580606742, 1.4951927148577118]),
array([1.706407852341273, 1.0322142875391824, 0.750398239048163, 0.4552658487623113, 0.0]),
array([2.5120515709427598, 1.7976413530107545, 0.0]),
array([3.136886169760116, 3.0056589617795826, 2.0854206541152696, 1.8555474459682344, 2.1045496710753135, 2.2367317774474254, 1.7348035962412411, 0.8035740200189553, 0.0]),
array([1.2168043523634675, 0.21318619266649574, 0.0]),
array([2.6117909237379164]),
array([0.7972629465657028, 0.0]),
array([2.3678783550812157, 2.310876131004744, 1.839988118951866, 1.435902118430155, 0.8353939962129314, 0.7811874853847285, 0.1702894774956929, 0.19904635048716124, 0.29526872429405826, 0.6710533009649298, 0.7423971475598997, 0.0]),
array([1.9311286121929885, 2.4267054818741327, 0.9031475705970772, 0.03374454515213719, 0.0]),
array([2.406693329896115]),
array([2.468070742316865, 0.0]),
array([1.047595697937218, 0.4628540906404635, 0.0]),
array([1.3869589682793162, 1.5570326464917095, 0.0]),
array([2.1730786833219877, 0.0]),
array([1.9963807519103653, 0.0]),
array([1.7957232606576168, 0.0]),
array([1.0494397682044703, 0.0]),
array([1.1846594422605887, 0.2797434393352335, 0.0]),
array([1.233052858755049]),
array([1.5286826857036198, 1.0572713762569197, 0.37065609687281764]),
array([1.2971877055735028, 1.1773719750346632, 0.5530868512615483, 0.455676191196609, 0.0]),
array([1.2852655155645805, 0.0]),
array([1.2065838223354353, 1.1377097277526307, 0.2275154858876166, 0.0]),
array([1.2025439928914285]),
array([0.8055449067313577, 1.2688233299024234]),
array([1.1575017053023422, 0.0]),
array([0.9820498371186495]),
array([0.2410136944207142, 0.0]),
array([0.5214323663396999, 0.6018369919751824, 0.0]),
array([0.39279967868905624, 0.0]),
array([0.2276636528165858, 0.0]),
array([0.22080305525826724, 0.0]),
array([0.02917790641248321, 0.0])
]
d = [data_1]
names = ["53"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T23', 'T24', 'T25', 'T26', 'T27', 'T29', 'T30', 'T31', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T42', 'T43', 'T44', 'T45', 'T46', 'T47', 'T48', 'T51', 'T52', 'T54', 'T55', 'T56', 'T57', 'T58', 'T60', 'T61', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T85', 'T86', 'T87', 'T88', 'T91', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T100', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T109', 'T110', 'T111', 'T112', 'T113', 'T114', 'T115', 'T116', 'T119', 'T120', 'T121', 'T123', 'T125', 'T126', 'T128', 'T129', 'T131', 'T133', 'T134', 'T135', 'T136', 'T137', 'T138', 'T139', 'T140', 'T141', 'T144', 'T146', 'T147', 'T148', 'T149', 'T151', 'T152', 'T153', 'T155', 'T156', 'T159', 'T160', 'T161', 'T163', 'T164', 'T166', 'T167', 'T170', 'T171', 'T172', 'T173', 'T174', 'T176', 'T177', 'T179', 'T180', 'T181', 'T182', 'T183', 'T186', 'T187', 'T188', 'T189', 'T190', 'T194', 'T196', 'T202', 'T204', 'T210']
def get_taxa_names(): return taxa_names