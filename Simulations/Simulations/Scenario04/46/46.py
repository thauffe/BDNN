#!/usr/bin/env python
from numpy import *
data_1 = [
array([28.47213874520204, 32.03983257986424, 32.27169760258688, 33.46855985189782, 31.48179174444117, 34.840472312552436, 29.08074001709129, 29.09094561331132, 34.060907179837336, 33.95935691023657, 28.875606347935417, 32.17657786939016, 32.718353485559625, 34.59501459020855, 24.017198629682262, 27.009214860862993, 27.11876073771867, 24.091507415841665, 20.618433628879814, 20.49780954678745, 22.64635611548298, 16.369779171852684, 17.569873927091553, 19.788479864584055, 17.762121327458658, 16.1821644970531, 18.91730015052715, 17.878030803327185, 19.752938955853516]),
array([29.70942319186758, 29.741619471832447, 29.096982231243505, 29.392980609746594]),
array([29.237463904434442]),
array([23.525057951399727, 23.747717699299006, 25.12099493490808, 26.25877468380419, 25.731634003114507, 25.54839837021222, 26.270068073326136, 27.061881247765744, 23.507648903754127]),
array([27.321119675710978, 26.218000116080816, 24.974644550868796, 26.67303359942287, 27.585316348761186, 26.271530283700923, 24.94647410500181]),
array([26.804945096660216, 23.77135981083447, 24.857604178809655, 24.73674077104922, 22.548500448055183, 22.504725110424218]),
array([24.677466703403866, 24.92660767689915, 25.834837755159235, 24.268190766251898, 24.2819385858498, 25.96584671775544, 25.694215669086, 23.04026189767658, 23.281390268595054, 25.403178468811014, 25.04965167991784]),
array([25.905684818949137]),
array([25.41841603787565, 25.635448290281, 25.518998689417852]),
array([26.19354415234553]),
array([24.058910916043477, 24.876611878223578, 24.29692907398959, 24.44372045726308]),
array([24.15471623505439, 23.395965569814347, 25.477913235231423, 18.749158917692125, 16.354231397469626, 18.931616785724493, 17.55844281028606, 16.62961476890266, 19.583895038120513, 19.558719546330803, 18.363446678117104, 19.304701900399294, 14.657230476970295, 15.350468649055431, 12.080232969658418, 13.006853936367971, 12.95486797802536, 13.732255725627999, 12.268521386248338, 11.712440000606495, 8.729369273653413]),
array([23.05046305787056, 24.4660988944841, 24.95178566262685, 23.96771918707245, 24.687116596600447, 23.831123009544147, 25.227163600846986, 23.21618324182328]),
array([24.4242119798437, 23.200456755330535, 23.37793261776641, 24.41488204861982, 21.039744150499125, 20.774355000512553, 22.340149774136613, 19.289865791281244, 19.6703589759797, 19.897939987239837, 16.304576801293578, 20.0760378231982, 18.42009262266808, 17.566563722161334, 16.8126878889164, 19.88921139561495, 18.893038140367437, 18.429887817053327, 20.365872050332943, 15.910793186127774]),
array([24.171403555436935, 23.805922467899048, 23.194771555825554, 24.581376183040387, 24.175839695475418, 21.35506759605072, 21.89720409790453, 21.006878822701978, 20.40292552591872, 19.742321760157715]),
array([24.634176987520977, 19.366207562381547]),
array([24.318858585681195]),
array([24.09513965974105, 23.38436774661606, 23.155947203239823, 23.91657874968353, 22.651013286607444, 19.82712464476661, 18.634809692979072, 18.881655240235325, 16.575866759354895, 20.418810499430265, 19.283934469254227, 20.405047519276945, 19.53729810344063, 17.820825607001925, 18.999747946451993, 19.475592149639212, 17.656393119478633, 19.35214393740121, 14.80787476802397, 15.561759748738488, 13.879760039920285, 15.324490342924038, 13.062189595769148]),
array([23.534579981352003, 23.971149699757547, 23.9152918448017, 23.496919822979386, 23.587188974233296, 21.11815723996506, 22.83085219001837, 22.912741888055923, 19.882720092636337, 19.63533193780896, 19.68170469534916, 19.855617445095852, 20.141773348232054, 19.90471042539341, 20.198525636639573, 19.80463799271691, 20.012092130655653, 20.41229582715412]),
array([21.660472335757078]),
array([20.73773028949801, 22.314468357587458, 16.989963918177725, 17.319670246851906, 16.16630073819154, 17.915955521963262, 18.8968927904394, 20.024778843894666, 18.09451600983597, 19.888896637132547]),
array([23.377693959603338, 23.425333783076248, 23.126451902355498, 20.716161371644016, 20.85925629396068, 20.66992872547384, 20.494918809402172, 20.54403282974868]),
array([23.080588888357614]),
array([21.54747402455121, 16.334309420537156, 19.262549627807534, 19.67503856555219, 16.896463769224034, 19.129285128361108, 16.583650856522215, 19.22378771540762, 14.148581370719352, 14.201717853670926, 15.281825307407162, 13.510446460714924, 12.330585873684544, 12.317127207554133, 13.008882805578915, 12.777563618524479, 12.088936127194028, 13.733226916094724, 9.199504467220802, 9.711963753134635, 6.153875977684377, 5.144281502821643]),
array([17.648050567298657, 18.997561640778763, 17.60666230356474, 18.828833784313822, 16.527990531480935, 16.259335488023225, 17.069595790455253]),
array([22.8249563663707, 22.652701147514488]),
array([20.51438620026311, 20.991535644176057, 21.349823059446308, 21.02417861737921, 16.93750579524006, 18.042125383289196, 17.846891315685625, 19.181944500925674, 19.103436971646854, 20.078175276459653, 16.479980832301905, 16.482986400652205, 16.37889402962864, 18.040082588824248, 19.848747438712522, 17.119027572624514, 19.729029951073674, 20.315519590169227, 19.40095074528129, 17.19802891116527, 19.696820145115563, 18.950078352267848, 18.12951638792073]),
array([21.58049472199518, 19.809518090343015, 20.37318575907655, 19.35529642979978, 19.964984434871845, 19.574185349359336]),
array([20.95473270345193, 21.814705656354423, 22.458527211291013, 22.21244113750662]),
array([20.55184265632809, 19.7559129665838, 17.621292094380454, 19.246903500672417, 15.379677362984056, 13.774836670840132, 6.836286408704072]),
array([21.566913346455237]),
array([17.00236215695116, 16.547899886785917, 18.344834677373083, 20.095062996774473, 17.25521780057578, 17.998220215817383, 18.610416391175256, 19.58731156321586, 19.454328887544207, 14.826630695507143]),
array([20.788028524664306, 18.784773668056232, 18.44356702048839, 16.81610678506517, 19.05006967502247, 18.53619051079317, 16.219226216424545, 18.792596601793043, 16.776659722764805, 16.397296298083695, 17.123667854163966, 19.08008643723938, 16.058666271394543, 19.9450982249505, 20.257892123928254, 19.145497120383855, 18.337955180082368, 20.0312218494569, 12.592449746917675, 13.123777691534261, 12.235319716569617, 12.76410948540088, 12.720249070182767, 12.92549293102846, 11.65430289646156, 12.535162904934989, 12.248002015793663, 13.166457084009435, 11.815026171795907, 12.623174667249302, 9.694005553445022, 8.83932094717302, 9.932669108590288, 10.696570428248565, 6.053427481422353, 5.467605261700131, 5.391166124737559, 5.87898181451272, 5.942332138035856, 7.058020409646696, 6.712841057500459, 3.7545703647969693, 4.262256755889083, 5.253733762961504, 4.011520224419518, 4.717298301976835, 4.289843149235815, 4.998849106340316, 3.205837192318073, 2.119835397618819, 2.4502885851680847, 1.6797738163708444]),
array([21.135884326341206, 20.614187177962215, 21.12925760848326, 21.405700250516723, 17.743381997002885, 17.743898647691818, 19.148666875816126, 17.35017414330391, 19.91448285477616, 18.07775243220893, 20.15394092428234, 19.461353039614348, 17.69720130190119, 18.330708836855084, 17.55736640487064, 16.721896550457267, 16.856360370203323, 15.944542175990248, 12.509954400936644, 12.696118128395955, 12.418803784394626, 12.226671641568494, 13.34226127266208, 13.70787982890447, 11.888407363832444, 13.742603112841513, 11.849379045792345, 13.478826731055713, 12.118592738680185, 13.717968620508632, 7.441540880756521, 10.942714490028365, 6.987942844871209, 5.811640321930941, 5.990729029295826, 6.972959815076781, 6.128152050615805, 6.814621646300753, 5.76636935743836, 5.778409319414218]),
array([20.760264467786143, 19.624073175276198, 18.27444436822076, 19.350512188406018, 17.248170132185784]),
array([19.858618748206585, 17.15360185060681, 19.47288271125827, 19.634887150454222, 18.942132130025776]),
array([17.70597552955177, 19.179331597956672, 18.341405225961804, 16.78720028228288, 18.44218453526145, 19.061372228281268]),
array([16.95801179596485, 16.975793153393948, 18.902274802650894, 16.582781490560222, 16.73158992918723, 17.169007244375045, 19.394042277690563, 17.067728535164772, 16.291835537420766, 18.641875743170342]),
array([20.62131429931683, 20.53502962928818, 18.86068163273788, 19.586154849687105, 18.798075083407408, 19.44157079764568, 18.800223044821617, 19.269229810071863, 18.919348885566794, 19.884165111716396, 19.085393358834793, 19.927912088365456, 18.508958191519746, 19.289620081962084, 19.026924465944372, 19.538243246701686, 19.83286883745315, 20.1502791307018, 18.70477530952249, 19.550744208607057, 19.21961036346327, 19.572851886467785, 20.042410612996417]),
array([19.708258666556684, 19.40407104864556, 19.325341743314713]),
array([19.342759648262813, 18.35633953637037, 16.082393401591236, 17.951568079736052, 18.472734698168175, 17.87126745063949, 18.29110899328695, 18.031471540921387, 18.45095661430131, 16.629918864498272, 18.674672639721642, 17.849417734737663, 18.867627514911124]),
array([17.863510679566584, 16.070301331908293, 16.267628281670984, 15.985164619432823, 19.42127490329297, 17.54774384958747, 18.806720206920183, 17.46517855618753, 14.155671339037175, 13.785644185769037, 12.054653812799664, 12.993639014292445, 13.487612244199063, 12.872197951213348, 13.41683071108708, 13.277217144094879, 12.667175637941053, 11.957866628970985, 12.68171766977063]),
array([16.514302260687064, 16.332844014827764, 16.825084007101196, 17.215535488015398, 17.838254524226038, 16.237229732219635, 18.529961913063104, 16.261799248772498, 17.391847619822244, 18.942825706908422, 18.700950869342815, 17.70248013338704, 17.02721284876352, 19.60575188580687, 18.058758385790707, 16.78292743743835, 16.335805328779138, 16.15727891121233, 17.635102254166757, 18.56077666954813, 18.741845170888407, 16.659535310793217, 16.920702172775684, 13.906365068328123, 14.95915670807107, 13.988302637546035, 11.908453047168123, 13.670944763030723, 12.787034532474244, 12.008530537660361, 12.181509104108914, 11.835581015410334, 11.794031461033263, 12.002454545198148, 11.898045259038463, 13.723807069545707, 13.721446148264462, 11.822718253018598, 13.318593093522063, 13.24935849478226, 9.429589141016352, 8.142895102135466, 9.203718082449726, 9.574181439035423, 8.725543774862462, 11.604698005299953, 5.548151940996291, 5.465151563551323, 6.229820100673161, 6.347391996934966, 5.956208593399551, 6.789985310306638, 5.485464382188843, 4.943573438650633, 4.554725574896621, 5.288222048032297, 3.992852451603319, 3.91629687813528, 3.743444520166414, 4.605555825772173, 4.596472687652087, 5.3116270233867136, 4.444967306800843, 4.920725122564384, 3.5659184974766487, 3.3312158962432923, 2.503592147580326, 2.3456065968503106, 0.5135697612359117, 0.03319850643394437, 0.11637925630975263, 0.030007657960018003, 0.05036656385376066, 0.0]),
array([17.621841497918854, 19.291253732186764, 16.52006819052284, 19.42009187502867, 17.89534089401857, 16.493941515166725, 18.870445492536877, 19.69856199757922, 15.515831590842183, 14.838038046293537, 15.218357467750163, 15.92690704264096, 14.62071735893309, 14.21211223824767, 12.484266294262321, 12.812582363238459, 11.69070484915419, 12.959848011435184, 12.626761939319692, 11.73323180829605, 7.615228800906946, 8.71935762382511, 9.916653562454972, 6.296870651793333, 7.053294618513554]),
array([18.538707222430094, 18.817645044895176, 19.07427134880182, 19.331620574542416, 19.084879853076288, 19.345596087211543, 19.113347124809074, 18.84465199850639, 18.782637718020354, 17.998579642282248, 18.51289361494924, 18.337132143308082, 19.036217183649534, 18.917138007047438]),
array([18.118349471359732, 17.479205557356472, 16.442859169119874, 15.990241434511642, 16.5400534970966, 18.000603591271638, 18.14803504067644, 15.607081887655873, 15.864858119951842, 12.058089347767298, 12.721686747813004, 12.274127851220964, 11.71260470324382, 12.920184975800943, 11.705760508432963, 9.742068750656045, 8.006429957344848, 6.961461874931492]),
array([18.043830226146177, 18.176558140151673, 17.567685544478106]),
array([17.01954462875893, 18.32273361506855, 16.87179511459622, 17.32565859756858, 15.976108414407495, 17.121846477428296, 17.51529389965414, 18.071204601801952]),
array([17.321455070524497, 16.677802709275145, 16.23191785399573, 15.599517191946513, 15.879323643712747, 14.65523936074212, 13.605285548234285, 13.480989512522564, 13.429546563388179]),
array([16.92919601359254, 16.47183888378262, 17.358301758669448, 17.543135026370585, 17.42266980741553, 16.570512967061998, 16.101056748793653, 16.897786652287856, 17.263408244956906, 14.44283072176333, 14.480521152885759, 14.394179067869064, 14.24428540862895, 15.591127921194897, 13.950199575329531, 13.280752548759157, 13.432422095530116]),
array([16.57186539767378, 16.918421790879396, 16.840572508660813, 17.513197188576665, 16.655749387436174, 16.460014852120246, 17.04415209191414]),
array([17.324048185965882, 17.26507700886054, 17.243756439059343]),
array([16.061866911554755, 16.25238373678144, 16.72875915233719, 16.28703872579845, 14.137735277037246, 13.822520757153177, 15.758665405426092, 14.829204612338083, 12.107520740083508, 12.252549549613297, 12.220042395492861, 12.65764130303775, 13.030763298022423, 13.096089873711826, 12.112205037057757, 12.680712103068329, 12.041129940648688, 13.704287839778011, 10.351908106517296, 8.308190876747755, 11.026131023495996, 11.22801512871565, 10.857956863656435, 11.257685397434436, 7.565660400997586]),
array([16.605828442814353, 16.35576596714564, 15.774203394321054, 14.374673669154596, 14.626127575922348, 13.906784088196652, 13.658850036650458, 12.782420281498585, 11.826687544289697, 13.07235930452149, 12.134480672106491, 10.658188476560394]),
array([16.668778477203645, 15.476288071338345, 14.90727899285447, 11.665867079292397, 12.209468871875838, 12.426014273446023, 12.673967204136993, 12.415159489837473, 12.420734702965559, 12.5346867271369, 12.786015793120091, 13.735335028342572, 12.249990565233373, 12.308552548957913, 11.724816294189399, 13.64111108612042, 12.214034870929419, 13.069905232239286, 12.631433850261999, 12.731176321825167, 12.880906390748025, 13.158538830753715, 11.56156043447583, 8.251210283091124, 8.83488910605031, 6.482016312997487, 7.200353590253149]),
array([16.07384067323741, 16.385348403620945, 14.75240280393168, 14.548635424325516, 15.952244822672371, 12.169626317495915, 13.301860532714528, 13.154915662159732, 13.305786836749299, 12.910358342260407, 13.061258419194164, 11.68897744317351, 12.70938205218735, 11.790117162378795, 11.788273100225943, 8.80061604901379, 7.586421252145181, 11.15669275166077, 9.132579061682373, 6.790833717211778]),
array([16.59016344611931, 16.038399770016163, 16.754406484854396, 16.40481819111272, 16.300956826186063, 16.608708135389396, 16.435825391356484, 16.08113260465555, 16.252954859571172, 15.124521282660282, 15.2169546808103, 14.17379449552056, 15.450655927275914, 13.20298082047511, 12.371481199308292, 12.759183414961596, 11.720241753935369, 11.89957615968482, 11.922245106864224, 13.677009148172706, 12.069144848624767, 13.45516629993285, 13.750367871556561, 11.73779723521027, 11.988268609885804, 12.407971872903644, 12.539914272059319, 13.21162069058611, 13.108446853668287, 12.829647779014786, 12.240692955509447, 6.904668965670099, 7.129266747832983, 5.000104884066484, 4.051344670304005, 5.003456249408396, 3.8190642972704056, 4.05548669145979, 4.88953905454423, 4.7341577918339395, 4.41085155027264, 4.475131879635158, 5.226014161129676, 4.08133508653482, 4.446790446625166, 4.708388646026474, 2.895964759025933, 3.516988006307924, 1.950904488931617, 2.195261836369721, 2.3428844196124885, 1.852875777992645, 1.9538474556802718, 2.4268418574843436, 1.472461538878141, 1.315947336775377, 0.41677242372348045, 0.7421960285559691, 0.0]),
array([16.525937720081473, 16.215819578696554, 13.904186365272597, 14.283731474078415, 13.83721144832135, 14.945908644195987, 14.073651505395267, 13.340558613261173, 13.413758016530759, 13.048932973553994, 13.037269644075007, 13.806749905777243, 12.621212061245888, 12.16536454665282, 10.874788413603056, 8.28712860393718, 11.431084922275016, 9.727888029839646, 5.426167389573161, 6.806055963852913, 6.898892872714995, 4.190946082013664, 4.5654556394251795, 4.30934580431251, 4.939111888571917, 4.394722101846027, 3.6168544445812802, 3.8145228737413213, 4.872122608387715, 5.073494231827411, 4.949228662381942, 3.8163436596571345, 3.3501874462790795]),
array([15.110894630844115, 13.805719651270824, 13.710305361579664]),
array([16.388138560389073, 14.92851224818139, 13.613199849606113, 11.879352450076414, 12.43472310154834, 12.156144595433485, 13.523442137055316, 11.649567793763348, 12.445624229532633, 13.076115932364978, 11.761392959668257, 8.487567175855263, 7.828516725981028, 6.79944187787047, 6.52406384489208, 5.209532488198194]),
array([16.249497836036937, 16.27243676919575, 14.923692641956798, 13.075203729510186, 13.09654712906593, 13.479702011255641, 13.356546100617287, 12.868437278867695, 13.380934060761039, 13.043871222220528, 11.972064063711171, 12.46149293270789]),
array([14.136905172627127, 13.395455389504956, 12.18272205976812, 13.10471605876263, 13.723126302860214, 12.66609912956932, 12.346708310715877, 13.438733520766855, 13.271851069091063, 6.320246151375886, 5.884926062880898, 5.862010429041391, 5.153566928641722, 4.258485344899832, 4.812367943301486, 5.10775382296912, 4.574096979384663, 3.817274882473215, 5.26912002820612, 4.0718854767306665, 3.5915877911379677, 1.4365773382378753, 0.3691767743324134]),
array([14.05591161755172, 12.350048975426697, 12.499212077118813, 13.582547475076664, 7.365317219211517, 7.052768421914794, 4.575488348460147, 5.27790359566446]),
array([14.958956401361764, 14.79127360531624, 14.413501024246148, 14.199737974570565]),
array([13.948084834869494, 13.177676122024655, 13.174461046316397, 12.538949096654909, 13.05541195090619, 12.285684749707231, 13.418555660891794, 13.031422374735046]),
array([14.8886726910235, 14.270093638693282, 11.902196361467436, 13.470839159657682, 13.255473509301096, 13.667862546648196, 12.056013291253326, 13.79402046220139, 13.754594804094465, 13.653923507387695, 11.832042999586653, 13.721376179031175, 12.338069323925332, 12.899879935411647, 12.677084827763647, 13.61187036600976, 13.775429505926242, 13.670668134225972, 11.490330175265976]),
array([14.908790524168744, 14.99229604606443, 13.394305981108557, 12.499411013168453, 12.918378515792474, 12.817748155032788, 12.275689561672737, 13.701099743680523, 13.11471499098353, 11.731237928850458, 8.952581724557268, 6.840702984072268, 5.149634983561075, 3.615657046592575, 4.523142126030491, 5.180424465248609, 4.07949105070321, 3.815360781191643, 2.250153261221066, 0.6547238244157995, 0.0]),
array([12.596404950342865, 13.157010203924765, 11.674117160949129, 12.996265459967294, 13.068020615282716, 8.211724174508133, 7.33802998140523, 6.427031498795499, 4.90166104282104, 5.1231139294974914]),
array([14.928117351452903, 14.232864897782388, 13.639928346281119]),
array([13.910143584236558, 14.492580090675869, 13.502999375263423, 12.342086987424452, 12.370066930074831, 12.187831606210363, 13.577479037381549, 12.854896506429473, 11.917518266785748, 13.405321746163747, 12.341195288918419, 11.642842289229032, 13.774787185660449, 11.975342135202634, 13.457282307191631, 11.878446222527364, 13.358735930801902, 11.746653192753529, 13.153329174082419, 12.128047238669105, 12.457875186334727, 13.224763814384275, 12.234678610848519, 12.713469820638657, 13.219330075925367, 12.082717459680484, 6.899934658332546, 7.21535401665397, 6.074718642886398, 7.229542316821675, 4.413202631773879, 3.833110434589493, 4.166079689584464, 5.2940104511104, 3.9804165794224486, 3.887169919383603, 5.09630851812627, 4.917671946241003, 4.621909508668247, 4.796500895802662, 4.2935464227025735]),
array([11.950970649703265, 13.261321027376509, 11.635559491321642, 13.20129004959728, 11.713495831415974, 12.177637227438385, 12.681924388393279, 4.898707489756735, 5.233573180881199, 4.808771603034998, 5.1554045718619745, 4.812652940286242]),
array([13.893419609575709, 13.241679971715055, 12.176766401048626, 13.313501173735398, 12.863309966919088, 13.695147623509465, 11.824014241605031, 12.00851018631718, 12.391966765828098, 13.050366727925807, 13.461046122211151, 13.45447258095503, 10.119496360405293, 8.368419434509306, 10.724293915181248, 6.0237376088057015, 4.8164544719127065, 4.388283638172071, 5.15937603888585, 4.46556625444442, 4.8939948052421895, 5.140385470040439, 4.4797156124994935]),
array([11.634976353380845, 11.849029735089996, 12.191314488113816, 11.842514223030255, 13.25484561126308, 11.901397082407598, 13.43079202302598, 11.827184897533957, 13.80598295977533]),
array([13.35132462935244, 12.97151869991122, 12.968855179147326, 13.712238989656335, 12.613199780600498, 13.173222044557662, 13.205022447599061, 12.957475288174344, 13.20644604692011]),
array([12.220783792913167, 11.995442920899997, 11.799192263950639, 12.337265433456006, 13.254131347368718, 12.535012847335459, 11.71765159294803, 11.918473748230408, 9.790346313973366, 10.696178856582591, 10.496583089076937]),
array([13.625338391552974, 13.527248620979218, 12.919623971035644, 13.58236724440354, 13.695553604832426, 13.36459074090287, 13.696896215887756]),
array([13.041723781250152, 12.309807977833497, 13.36912930145727, 12.744851406719759, 12.266550534953737, 11.852405502585771, 12.866476385573973, 13.508447165875678, 12.42337594293937, 12.26154279271253, 6.2848267608421615, 5.972320269353835, 5.198299191921115, 4.004110391446692, 3.8010023661744965, 4.702521142384169, 4.446451388066373, 4.9249431696553865, 4.208267713850187, 5.214624827162628, 4.124883733005628, 3.677235277950416, 3.7104229825431765, 4.356504095023211, 3.0616520689908286, 2.7798553720601804, 2.4027098693352236]),
array([12.959395297584052, 13.450358697762697]),
array([11.712375311389378, 12.962373885718343, 12.636757200176023, 12.211778787669447, 12.063965014217118, 11.662862578470502, 12.360807576499813, 12.223132880851828, 12.328374494744248, 12.732363888560906, 12.529300560363863, 11.881586984807369, 7.700068709464643, 7.946679984865288, 10.489794253529528]),
array([12.615845617284547, 12.758247959520906, 12.814661161079973, 12.379999379581015, 11.956362517792105, 12.547242046835736, 12.535215842106854, 11.892349065453462, 12.827166905172826, 11.659742423848623, 11.962726283150982, 11.854171236735748, 11.967035663892709, 12.176995195703755, 11.897706241052566, 12.46633698127312, 11.20683829747435, 10.117693345134736]),
array([12.014891488060506, 11.916687226434624, 12.096663313937258, 12.396962427341707, 12.604052477652639, 11.987634393739166]),
array([12.623813254224237, 12.534145154491224, 11.912726257952727, 12.95495366237579]),
array([12.462783133210138, 12.199058283669318, 12.848980251377418, 12.32246276847457]),
array([12.601246097400347, 12.433287452350893]),
array([11.732682477700001, 11.771262181398997, 7.737950178434479, 5.603865830724327, 6.292713447792441, 3.7892721085764514, 5.071686827708677, 4.051290097972499, 5.120224707773475, 4.525788709437373, 4.117419523103044, 4.124787570251208, 3.8119414795974844, 4.9710751391943635, 4.842754621015554, 3.889063571131584, 3.982207972812588, 4.176142942124163, 4.5783520937055755, 2.3468266530893955, 0.9154130060346135, 0.5501609546623389, 0.0]),
array([12.30674212366117, 12.01226090303265, 11.967681363764038, 11.955136955255025, 11.215429102365107, 9.626917333775904, 8.808840266829575, 11.127669308366958, 8.653473362602726, 10.35507732530551, 9.641802034936184]),
array([12.245013543817107, 11.859234522845897, 11.6658696990889]),
array([11.73083391928844]),
array([9.639116709141787]),
array([11.040029040542569]),
array([10.3588127752202]),
array([10.450075687475014]),
array([7.258284428691963, 10.400302477997956, 7.138857061250435, 4.147977682067865, 4.313548926698066, 3.6736600000169504, 4.972797131634708, 4.996655122265743, 3.7730711540591244, 4.741141500503015, 4.7783331262639255, 4.0376785451260995, 4.259472677366637, 4.934376915518158, 5.295921710429303, 2.6888519836016758, 2.1493844771613837, 2.5427148289308823, 2.0461856548919855, 0.2627820330783702, 0.04724660979330947, 0.0]),
array([7.925038516276034]),
array([9.637993296793908]),
array([9.972583862042141, 9.557572573786643, 9.709997050913165, 7.251402516297555, 5.043315606214709]),
array([8.908141360352095, 3.4718473476156406, 3.1180070387782153]),
array([8.9299659530217, 9.446910371855816]),
array([5.868020096458722, 5.257141069516311, 4.820268982760238, 4.459884930968597, 3.9480649493417443]),
array([7.060017854674525, 5.8096652298673375, 6.085693691512606, 6.849712472561131, 5.978754421539763, 5.3171341074387675, 3.663226763014185, 4.347231362915917, 3.9148029521881242, 3.929765640180107, 4.138689834322454, 3.7813158506441633, 4.136165660925537, 2.3375749312140335, 2.1488866982163497]),
array([8.825353666885743]),
array([9.096748346865526]),
array([5.662203885842814, 4.6150459716158325, 4.012084320032291, 4.14165795339442, 2.794610587386698, 1.9786769526896932, 0.0]),
array([8.55798660493164]),
array([7.815831290380922, 6.8887744679252005, 5.739676042680829, 4.950682826469562, 4.734021806953462, 5.119380254660633]),
array([7.077119280453878, 3.71155804037286, 4.45588741876439, 4.307509595181526, 5.052116834915913, 4.295396194881895, 4.532436112047071]),
array([8.057948285248647, 7.4596833517098835, 7.538365940470262, 6.260674432265008, 7.186443257021634, 7.172571828632288, 5.744345436152562, 6.970842686100751, 5.395302576785604, 6.40686664666288, 4.527129520461962, 5.025397600176039, 5.227971403085772, 4.598825717092262, 4.8845560134301325, 4.494205242849968, 5.24852363550387, 4.569889080057684, 4.272331072821984, 5.299273093263604, 4.609081504359045, 4.693662058698392, 5.161092424220641, 4.919159022886468, 4.46622677046226, 5.170031325932967]),
array([6.555489636445985, 7.080711307809324, 4.665813973492452, 4.546005351868509, 4.621636183969563, 4.829502702582649]),
array([3.769954336727504, 4.129044622337069, 3.981119283112193, 3.3548789636421623]),
array([7.016605532580308, 5.455841216492039, 6.14304389420418, 4.905758141510451, 5.020169394488214, 4.452451323284707]),
array([6.912397665345532, 5.382818693097904, 5.697904207690488, 3.760758305987122, 3.8605848311264257, 4.283566162898832, 2.3131605668304127, 2.551759510242058]),
array([6.13143732390971, 5.837536534126304, 4.424147544273963, 4.480736469619977, 3.728247658336651, 5.116272199090258, 3.994106698823812, 5.123722933834148, 4.510353821369817, 4.575346974225849, 5.145793573648475, 4.064211740445121, 5.051025465921876, 2.6389024719375476, 3.244493213733139, 3.4439091862365707, 2.9741205686515784, 2.1767130263296526, 1.810321142696392, 1.8010694896726251, 2.3213840583573155, 1.0946321601142146, 1.4442365443958525, 0.4535344443479644, 0.14241810650501274, 0.7772070036687674, 0.0]),
array([5.712498406167458]),
array([6.257967399286303, 5.818440229335214, 5.081392449451338]),
array([6.574558347918776, 4.415410507627614, 3.8214660736685477, 4.876537762784757, 3.638750254978328, 5.310536445598433, 4.546573608220901, 3.976277478967411, 4.7888660883008365]),
array([6.328662059580156, 5.351447270550832, 6.131924816971962, 4.8994583375475225, 3.8327588023932115, 5.042071091828152, 4.219471564437896, 3.987104006614609, 2.7123521627701637, 2.4324503015753756, 0.0]),
array([5.945728310910651, 6.014243865607778, 6.019468073341771, 4.87181566969081, 3.921683741826171, 4.977448871507306, 4.173410796082576, 4.0693401975321155, 4.725936673578276, 0.0]),
array([4.345842335505981]),
array([5.882270202446788, 5.439754694830185, 4.314367493247185, 4.540056790297062, 4.924549104228537, 4.21653732509813, 4.523582546209276, 4.884624941369375, 4.886052120038444, 4.4461298983000725, 5.030260160867427, 4.409247478685575, 4.317593313876346, 5.122613285305343, 4.764609401268261, 5.309165817966229, 4.167622649246713, 5.063145540355997, 4.769311876523376]),
array([4.749660143234064, 5.238742453111738, 4.37520491072875, 4.932253075029766, 5.141257031731952, 3.923000998838263, 5.253976605949952, 4.875499905118611, 4.590974837678323, 3.79750157971967, 5.2767280527969875, 3.4395051739783913, 3.5985673867140906, 2.566849048794616, 1.7137078330024083, 0.9651694131033214, 0.40348219482964887, 0.10855316310172361, 0.0]),
array([4.959172079676934, 4.582431429658149, 4.669031975444518, 4.876655650114488, 5.106901705377263, 4.851468558508842]),
array([5.686537438808445, 4.2422587509585785, 4.718659767890737, 4.231900598566317, 4.651603671522769, 4.187247079781439, 3.1674974038834263, 0.042858992558355805, 0.0]),
array([5.30126257608558]),
array([5.443524535902512, 4.948371273551429, 3.6568796364132163, 5.315448279479692, 4.577557466128057, 4.1563725961916, 3.6562301661919108, 3.9408606509210147, 4.952830499627458, 3.6193663906208435, 3.0890892077765635, 2.641791412008515, 2.1336993501809935, 2.455622238824844, 0.06927570738611535, 0.0]),
array([5.1574454623151595, 3.689637722280731, 4.320900354187855, 3.8011980216333088, 4.353277485952132, 5.200537367623187, 3.709660011477883, 4.071571235036893, 3.6527488278079088, 4.736336057089922]),
array([5.256914185441791, 4.184568185883403, 4.3152079878429, 4.328872687936343]),
array([3.684719993960809, 5.142719241382951, 4.317910228873069, 3.9312742363648177, 5.117216582200842, 4.555534928250241, 4.389462299560469, 3.880083680754354, 2.9119998318482105, 2.8181571576813926, 2.2749584755319923, 1.4417472972659136]),
array([4.476246371759532, 5.128238641251683, 4.60995418613343, 4.110252147653701, 5.033536527885375, 4.437275462883607, 4.1332726109589215, 4.54982294087936, 4.304525504279434, 4.7384907090171655]),
array([4.100889657974577, 4.686657413134356, 3.9562511101121816, 3.915072565988523, 4.90310118745063, 2.4543595088409735, 2.5691619203990053, 0.0]),
array([5.095277379896282, 4.570282746087704, 4.244712763775093, 4.3345829042383075, 4.976679952247681, 3.9257588106109846, 4.876082719976189, 3.818988017463135, 5.016396147783763, 2.3433599835080257, 2.3442938087245206, 0.0]),
array([4.614844103189333, 4.698203710675233, 4.154133122654487, 3.2554139727215112, 2.0459160652437283, 1.92114360404961, 1.040195503330144, 0.0]),
array([4.083773602077285, 4.368381365831578, 4.442585298022472, 3.928312456431954, 4.210088426379011, 4.685886827505344, 3.3319755241357676, 3.1491219662119123, 1.862549635847909]),
array([4.24101558330513, 4.5311159336543145, 4.4918375325787965, 2.984165480551369, 1.915614582315905, 1.466633720568066, 0.0]),
array([4.15546849953867, 3.864060423674001, 4.103103587453459, 3.6675043714581643, 4.161546465313332, 3.748384639725587, 2.843870645567591, 2.3031062910494864, 2.0189822637318477, 0.15644076661562167, 0.4350624019134154, 0.09683357478312687, 0.0]),
array([4.138974623722467, 4.02528931647347, 3.8549962254629757, 3.969629495793539, 4.366798669460487, 4.360028894132891, 3.8622654784925308, 3.541126268681579]),
array([4.313346462565111, 3.9527719267213373]),
array([4.1918865393149485, 3.7409289330165083, 3.9428536961112224, 4.2615843196173, 4.04453085719746, 3.7072317527801384, 3.6620108621338616, 3.689811453435551, 3.2658932931104614]),
array([2.6123788061105095, 3.243643574433531, 2.5239085221958675]),
array([3.8481078816282395, 3.6148853435181865]),
array([3.859743803937587, 3.4957846427721164, 2.9876602915542234]),
array([2.5934113380167707, 3.226580571429035, 2.4727794491707065]),
array([3.026322291030388, 2.9540729927341864, 2.3661755093576]),
array([3.7613292272921153, 3.281268166206719, 3.2713805262618223]),
array([3.226911093253508, 2.5680280857485522]),
array([3.3443363060008684, 1.8840007348972445, 0.1028429045767003, 0.0]),
array([2.9515178559160073]),
array([3.280688790088505]),
array([3.0406412266711538, 2.5262142946157855, 1.2458780365883193, 0.9948493561447854, 0.0]),
array([2.6373922692300362, 2.0561739047520526, 1.7790297833070357]),
array([2.180335178218017, 1.9671884941220819, 0.0]),
array([1.9430931296679121, 1.0577812686862829, 1.5550273022178471, 0.5262003053574564, 0.14476147436462805, 0.0]),
array([2.345394039772311]),
array([2.0401935920171064, 2.3950714301037364, 2.2763393455495162, 0.0]),
array([2.041218067775076, 2.279526133765541, 2.186772724553315, 0.4279310231023422, 0.6339286249167236, 0.29558779459446116, 0.07354112219420524, 0.0]),
array([2.2845937255857667, 1.9081063122402133]),
array([2.7060005148961443, 2.219986845498784, 2.016323422199448, 1.1185252469331541, 1.6110070888747856, 0.0]),
array([2.5547574070371333]),
array([1.8805210849614356]),
array([1.3420194840752306]),
array([1.6082276237153759, 1.1730578650282784, 0.5750327212068305, 0.0]),
array([1.2899941893612819]),
array([2.0981320903473892, 1.4567514549543772, 0.48425476043728916, 0.07517120327598736, 0.09094065653289403, 0.0]),
array([2.2051987136986493, 2.2234065477715625, 0.5711626587200532, 0.07039151276538874, 0.0]),
array([1.2509908667432774, 0.0]),
array([2.424782669137477]),
array([2.499401370983187, 1.967847024623712, 1.695974731807068, 1.478631545575647, 0.3440906111365216, 0.2572692853692533, 0.2153168106079778, 0.6244064517577518, 0.0]),
array([1.1608467844941974, 0.12564274782112037, 0.0]),
array([2.0134421502637427, 1.2963381183876796]),
array([2.1504248296824784, 0.7732269944249536, 0.0018883073722552185, 0.0]),
array([2.2183595017945534, 1.6539192901626774, 0.0]),
array([1.4462434485841704, 1.3108529236678064, 1.409434317112641, 1.2339290140687482, 1.3873073587841647, 0.30387611453885005, 0.2042296819825301, 0.0]),
array([0.25720208463602234, 0.0]),
array([0.08630207984256988, 0.0]),
array([0.2216073546746441, 0.316116052722551, 0.07379841870764439, 0.0]),
array([0.48130277830379165, 0.0]),
array([1.4160574784915954, 0.18124683609705616, 0.40696978259703326, 0.03656648585705087, 0.0]),
array([1.3770292125367984]),
array([1.4649543534204743, 1.297941541880915, 0.49322006802086665, 0.0]),
array([0.9287323035199383, 1.3106489580701948, 0.34615610994516555, 0.3410628374994287, 0.7520688828987819, 0.0]),
array([1.1081174090272556, 0.9158881061727517, 0.6588659734531039]),
array([0.7768556767459026, 0.0]),
array([1.0036196368745336, 1.289423569790722]),
array([0.9787586588498611, 0.6431375123138489, 0.0]),
array([0.8222964773176611, 1.297542324096544, 0.1089876592073303, 0.0]),
array([0.37143651452198934, 0.5904745170533261, 0.0]),
array([1.2263197555420329, 0.0]),
array([0.05325157770950148, 0.025629119942965475, 0.0]),
array([0.03522106463720326, 0.0]),
array([0.4823047855798363, 0.4099533318836615, 0.3904711194736739, 0.7709903130919689, 0.009409280302882739, 0.0]),
array([0.02643741155886098, 0.0]),
array([0.36578067979770346]),
array([0.6017027442842383, 0.0]),
array([0.5059274675788225, 0.0]),
array([0.2930085099703302, 0.0]),
array([0.4881150753123116, 0.0]),
array([0.2615392532346708, 0.00802318426628422, 0.10645652013797845, 0.0]),
array([0.4216514877514999, 0.0]),
array([0.38905504264287327, 0.034677739750888364, 0.0]),
array([0.09101991405631923, 0.0]),
array([0.2060051793942848, 0.0])
]
d = [data_1]
names = ["46"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T5', 'T7', 'T10', 'T12', 'T13', 'T14', 'T16', 'T17', 'T18', 'T19', 'T20', 'T22', 'T23', 'T24', 'T25', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T39', 'T41', 'T43', 'T44', 'T45', 'T46', 'T47', 'T50', 'T52', 'T53', 'T54', 'T55', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T62', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T86', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T93', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T100', 'T103', 'T104', 'T105', 'T106', 'T107', 'T109', 'T110', 'T112', 'T113', 'T114', 'T115', 'T117', 'T119', 'T120', 'T121', 'T122', 'T123', 'T125', 'T126', 'T127', 'T128', 'T129', 'T130', 'T132', 'T133', 'T134', 'T135', 'T136', 'T137', 'T138', 'T139', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T168', 'T169', 'T171', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T194', 'T195', 'T196', 'T197', 'T198', 'T199', 'T200', 'T201', 'T207', 'T208', 'T212', 'T213', 'T214', 'T215', 'T218', 'T220', 'T221', 'T222', 'T223', 'T227', 'T228', 'T234', 'T235', 'T236', 'T237', 'T239', 'T240', 'T242', 'T243', 'T245', 'T246', 'T247', 'T249', 'T251']
def get_taxa_names(): return taxa_names