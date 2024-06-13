#!/usr/bin/env python
from numpy import *
data_1 = [
array([31.774930182402304, 30.286961679858116, 31.740581105390827]),
array([34.21900445249469, 30.34762062448357, 33.02009249062976, 29.579159408801615, 31.282367180970308, 34.083564704298034, 28.141132760589784, 32.47261364219637, 33.14026682821967, 34.242223132425394, 29.261259570436483, 29.288234366469872, 34.20214919367883, 28.387138362944427, 30.953744064387763, 29.970997220932638, 29.927665391912925, 31.66507455611961, 31.38488787427205, 32.21399561611144, 33.827682326075994, 32.52357973761926, 24.187919626367155, 25.62760756790548, 23.755754776453998]),
array([31.107563721711397, 31.98998468657487, 30.313939012453506, 29.216009309604008, 25.922271580567525]),
array([31.46999626960753, 29.29052937613301, 30.562207963328387, 30.032467145831784, 31.602214826559358, 29.03615239172821, 28.70676257664236, 30.40839055176505, 29.17526170851304]),
array([29.498401568443395, 28.126124940921756, 29.172977574414542, 29.07362768481231, 29.303917156696162, 29.081691159780007, 28.142247779131285, 29.82278819935443, 31.450675654929356, 28.37987633062442, 28.209066964422284, 28.6076900702883, 28.695464086732727, 29.23714384003174, 28.441067546187078, 30.476396144321455, 31.40125180865939, 30.40158883983367, 30.48812520342169]),
array([28.551875173101166, 30.553246661492093, 30.577673109741088, 29.54389946736517, 29.135305246455694, 29.8771409399063, 31.129035471931164, 30.51676553666053, 29.686439654271194, 30.255739888440786, 31.22918591957836, 30.092867636093423, 30.852869295257683, 30.236377686194306, 31.340350772709545, 28.62803707695491, 30.163467380814996, 30.607076164536952, 30.737356128924038, 28.60446866488706, 29.083414543759403, 29.750830852693166, 30.30709679592331, 29.43732034789307, 29.988415562642484, 28.823414638330963, 29.606211718200786, 29.773846544559447, 29.131163259354686, 31.423525473574646, 29.993747755145012, 30.722238913784548, 29.286789689006756, 28.68411598468152, 30.807009496817134, 29.127755845151885, 29.57727368106733, 28.569916475261874, 29.978146039704765, 29.229368858638146, 29.0287193591856, 30.836680979700276, 29.346273972310975, 30.005818228950524, 28.71249810194997, 29.6744615813961, 29.7277831871865, 29.83538236766555, 30.453983464509555, 30.085626122548415, 31.188985790917666, 29.51160146371868, 29.774843376394973, 29.412898173188164]),
array([30.63687408084581, 28.68314148818523, 30.001838542842655, 27.191158334343484]),
array([28.327596419001548]),
array([28.80509466500712, 29.23481284727097, 28.787147530398492, 28.9380090779271, 29.154477968402894, 23.06592698628473, 24.30685794352634, 24.398505553916905, 24.678149364787245, 26.597657123791812, 20.619004489825663, 22.12714964656501, 22.11371990142917, 22.254086194336335, 22.54550486017849, 21.919518135283035, 21.173415881248538, 16.722359832549234, 18.09685498146623, 19.426155677135952, 14.692714525391416, 14.719746764444187, 15.462188116720665, 11.773260574473952, 11.116924156190306]),
array([29.03741434117448, 28.34159350209872, 28.124282844138918]),
array([29.152582002610362, 27.37593390033031, 24.006410250462775, 22.851863775783666, 22.855768580873193, 20.449238435550228, 22.65617312693164, 22.39923210655931, 20.656228284228362, 20.802090965293498, 21.35414116471784, 17.215279382409243, 19.598235507698057, 14.249767590555646, 14.793265716851437, 14.197063631509343, 14.020601770362203, 14.92498041197934, 15.641771028825806, 13.91077972646898, 12.919551567205778, 9.689065207925045, 11.116893332575188, 9.420031042769397, 10.065168856973706, 10.635533559698693, 7.8902491848002825, 5.8258869652875305, 4.668751021518762, 3.0465926702904813, 3.338423354179469, 3.544364527136403, 2.9165616146483333, 3.487576376918288, 2.6682458005583167, 2.689239607185199]),
array([22.51499346638382]),
array([24.06075089522316, 24.867134190521636, 24.280046329643696, 25.19855878001105, 25.359974904195337, 25.15359598032416, 24.922305419629414, 25.267177512291543]),
array([24.878611403988707, 25.041522074682366, 24.947131649319658]),
array([23.40940754319217, 21.953899990363766, 19.24948736089769, 14.871470352243076]),
array([22.505387185434536, 22.580261406908914, 21.098869012186167, 22.33795308139335, 21.006846616193013, 21.38735952931896, 21.148601804613985, 22.676407003361728]),
array([21.015014711077, 22.406166910456992, 22.426913002292256, 22.50707802174219, 22.44611547456883, 19.591741409451938, 19.83501444387798, 18.100012588406898, 17.60454118664367, 16.409549981989983, 17.210560591267125, 15.764063410605166]),
array([21.839982413041383, 22.32966524966422, 22.732509391940674, 22.030940882793203]),
array([20.930946990004347, 21.074300180790157, 20.50314905476996, 20.707348095136005, 21.043540302922256, 21.184959002387732, 21.31144897732324, 21.444826024352928, 21.323652329558207, 20.136585594340335]),
array([19.322783477675358, 19.86737834195142, 17.979831144168447, 17.50618605830607, 15.21890621284483, 14.126035058340427, 15.199676118702715, 15.915168319644504, 13.896721551822075, 14.249179919626021, 15.077449543801286, 12.4078061682123, 12.843934387647325, 13.64031594935558, 12.22409393644344, 12.878884841504615, 12.493807380572225, 11.161855270814623, 8.63644612483192, 10.911470974104258, 7.516236895511402, 9.1855762827363, 10.139752673679121, 9.447332855476834, 7.238466871221111, 7.045653361162654, 6.5073643196665305, 5.728756822518486, 5.400344846604256, 5.719550286566226, 5.996893162469493, 6.897197159819488, 5.028450746634577, 5.329003679113349, 3.4270457509489387, 2.7621509624251868, 2.68824067091297, 3.4189159564839433, 3.4118607862900907, 1.8382583226356142, 1.8461642864722396, 2.0935235714005476, 2.2560252578231017, 1.5396541359091638, 0.862006364790108, 1.3537333846521138, 0.0]),
array([19.15621602239774, 16.379094244787776, 18.436651290035105, 14.433489924200373, 14.35657979407129, 14.823862758974002]),
array([18.9379409246671]),
array([16.905984934168423]),
array([16.92705449738079, 16.20988148442583, 14.052729452591986, 15.140287700222954, 13.822403740148454, 15.704923726465386, 14.36395770649059, 14.297895665042912, 14.639695640690258, 13.926108942997194, 11.665732120964703, 12.215847198559988, 13.697653971398669, 12.79722761469812, 13.36509675239031, 10.707244249227973, 11.13274448242082, 10.250196772425923, 11.381862537487644, 10.013398125472944, 11.001695894998358, 11.258661889495633, 11.29175426722743, 10.236273426585335, 10.435605744366521]),
array([17.930333774472327, 17.291512812369504, 16.951255074705486, 17.851621011226953, 15.247231472374832, 15.577675275871252, 15.394385585751111, 13.678050785228379, 13.59484360879132, 13.789329858170229, 13.494136134566727]),
array([16.90137995285756, 16.837023448420016, 15.694866174074257, 13.849817890208913, 14.739647067590377, 13.985117232233863, 15.547359533231187]),
array([14.651707855010114, 14.132885582222142, 13.997321984744769, 14.316455166375722, 15.68001567461332, 15.03961537639378, 12.816307389307337, 12.42921362840574, 12.750706972193383, 13.664112228314263, 9.66808619239507, 11.11094705705273, 7.956483446966349, 10.907270190739407, 9.411604200286668]),
array([16.505693459079737, 16.341463706284987, 16.02079342781798, 16.45091011812565, 16.122884903070513, 14.15843466379024, 15.514336261170364, 15.404487188982888, 15.168849922003357, 14.190729978642066, 14.931075058863495, 13.98294419302307, 14.459372779830696, 14.795704419741499, 15.681610807907099, 14.500648831769803, 14.057784508941275, 14.842786817583713, 15.88089696224086, 15.560903322259866, 15.179833033154123, 12.355034841131037, 13.565188166168372, 12.28344094989205, 13.620655222961004, 12.575446314540988, 12.369397712211281, 13.206885685800106, 13.588019060749138, 13.387479511275178, 13.819863099067735, 13.776728429770081, 9.955814508099895, 9.217642800913564, 9.736112140227275, 9.674190308251552, 10.239949427771242, 7.523502893680425, 9.643261575034723, 10.549944275094072, 8.756913910138383, 6.1661225072665395, 6.5719048902263335, 5.981896114393143, 7.094100166327848, 5.5127971211359705, 6.120831785166123, 6.8192704185983875, 5.667991421608135, 6.864683488874991, 7.12626805972428, 6.9152652576508515, 6.480657690839234, 7.176331870274127, 4.7155009758585695, 4.096636116401742, 2.9049795253208424, 3.2813674560228803, 2.9144617605348535, 2.895471965962074, 2.603194226883086, 3.2059133042458208, 3.0975517373332404, 2.133242909041413, 1.804400941819035, 1.9028695653981504, 2.177567247935179, 2.0044711237699486, 1.7538490900110368]),
array([13.114140903650236, 11.244425162444633]),
array([14.801073276840029, 15.570182916526933, 14.80310225827628, 15.244245342062381, 11.978288884337683, 12.911445515166573, 13.745424377177951, 12.87285598872666, 13.259690426471762, 9.9495511244624, 9.92324925180949, 9.723690440377116, 9.875305912487356, 10.140976937307963]),
array([15.033467852841264, 15.513389655505874, 14.509917039438918, 14.50805482136628, 14.59614436876303, 15.094156320818758, 15.112393367312876, 14.725666496135327, 14.567149968361772, 15.126109396029205, 11.776107176001085, 12.883198649036578, 12.03473310816503, 13.007236834111144, 8.94130749195017, 11.314627039369418, 9.116070783744572, 7.4205730665766145, 9.780613072226444, 8.120543742266452, 10.947240567917586, 9.219694181325474, 11.547456127821492, 9.102348189877418, 6.372048477682027, 5.9270806935013125, 6.874248781791289, 5.970960840627168, 5.413220093090762, 6.356343393173799, 7.197818412385459, 5.726601438107055, 6.832949929267581, 5.2374079096636565, 4.129311114487573, 3.9277094594169832, 5.229253362859097, 3.4411624784061425, 3.118996235949167, 2.664811440592251, 3.3182197242659583, 3.5385900247992836, 3.481243416208332, 2.6069205553169432, 2.571276921030138, 2.214237264767701, 1.806951527304891, 1.5335527502489226, 1.2747016216861764, 0.20956518842812555, 0.24988001476320676, 0.02401320021105849, 0.0]),
array([14.362943498041625, 14.269445386513729, 14.370412706578941, 14.536931660867952, 15.083896742163422, 14.696890329382837, 15.297475620351712, 15.087252097471334, 14.383542046162715, 14.060884522581482, 14.833129402248444, 14.638364858524888, 13.88002288139726, 14.473608469893462, 14.208394287258399, 15.546749658154752, 15.25415819176167, 15.463686380585328, 15.20527053211067, 14.384545970392043, 14.712060660960386, 14.064876701507082, 12.23361183627317, 11.938795035590784, 13.718562362859627, 11.828647465862709, 11.986641532335218, 12.777848058716991, 12.331674570174165, 11.933891771257608, 12.622139377921192, 13.031722241329684, 9.99554131693617, 8.230800454459041, 11.147856567265967, 9.62913085404975, 8.91444777227682, 7.902872729857453, 11.396065091220065, 9.209969515437932, 7.515019638408955, 9.893912998205659, 7.297168299371512, 8.840364326508777, 11.57797533973558, 7.568205264612648, 10.493209363673566, 9.192412348053018, 10.687626285972978, 7.649643603380798, 9.047573521151689, 10.121425701262428, 9.184277704556028, 8.798744429466051, 8.667576329446161, 7.508239044186553, 8.202797710169733, 7.359121521522201, 11.142731530695396, 8.28167045705081, 6.700813956759159, 7.049422252120476, 5.83847336495484, 6.249850516697328, 6.767987689495027, 6.917177347695759, 6.629501710640286, 5.926572010661768, 6.432573431521305, 7.056648352486089, 6.2391985978148075, 6.489105893127286, 5.71524586856262, 6.286741737381719, 6.935564261531011, 6.288406807487276, 6.176550804384879, 5.3705107466573825, 6.363505740044783, 5.729760632426261, 6.867052521689493, 7.159741782906927, 5.495449677901248, 3.6412859723635504, 4.260982217674668, 4.792361385188448, 5.020130208834702, 3.9858484634789058, 4.905596997521343, 4.383573344745963, 4.214376169469972, 5.046787773143066, 4.1670414366431325, 4.935450905215543, 5.134866315794644, 3.5686759619121746, 3.4995919471205084, 3.446470454414898]),
array([14.847908930537367, 14.971159872242625, 12.733852260860223, 11.233911464138783]),
array([14.531911101855197, 12.480252619363219, 13.719197827006692, 13.65843270392838, 9.374538587346054]),
array([14.648032165549672, 8.50015518874799, 11.014108828873688, 6.522040463233778, 5.855852539953341, 5.394197638574361, 6.972362912427189, 6.677557015502645, 5.234308429962437, 4.268984036761427]),
array([12.735143254723331, 13.364814790017974, 13.69986347653736, 8.773602335521037, 8.339884518128004, 8.628485086264758, 7.7354947305588535, 8.934148064787939, 5.767160372830073, 6.995461703668841, 5.091181292692391, 4.468723793446176]),
array([14.260832914864976, 14.41166051975996, 13.859056695428817, 11.973157665467639, 12.500715644697268, 12.879186550572548, 12.159754796613862, 12.269981736093326, 13.533230097403692, 12.364692268599738, 13.189524778849314, 12.183630551186848, 10.064955610433255, 10.199179084849554, 9.079998432756788, 10.13929369125559, 10.798884235157178, 9.642681691805937, 11.625253835120352, 9.922654640886494, 10.787686781535655, 10.412589367219763, 8.370259591092916, 9.436182631145412]),
array([13.55206325049634]),
array([14.09130674637186]),
array([12.393242858448833, 13.08649510043261, 12.977259447746542, 10.356807560708694, 8.765025311251538]),
array([13.065413428701506, 10.260059021645992, 7.136594998557084, 2.951278286777417]),
array([12.778193881652694, 8.57293081895587, 11.078826050282094, 9.591986689405335]),
array([11.6891642328893, 7.042158407176933, 6.540938230470555, 3.0470837167435767, 3.2475540128368108, 0.0]),
array([12.594678008870092, 12.206634428200005, 9.339529164261217, 8.257421351030676, 8.554944244444444, 8.868819992363687, 8.072088699237003, 7.6755525368243465, 9.299315651482626, 8.410433869049251, 10.478677460048962, 9.67647252189917, 9.315082709113335, 7.140080613496752, 5.347372140295025, 6.01535632758205, 5.443656790261033, 5.8416595062424514, 6.892395672347332, 6.966113071813618, 4.090646171617341, 2.7034496055017145, 3.2135632078780496, 3.3188342189263627, 3.3644302790758034, 2.863958852592614, 2.7337069281245467, 2.3341012826126564, 2.575052423034404, 0.8731557569738664, 0.9679007449314067, 0.2772338906284858, 0.5219870040358474, 0.35407248621260384, 0.08839514259211992, 0.0]),
array([8.980508295644952, 9.290647743621863, 10.703332381137681, 9.239483481906614, 9.811010997200171, 6.694199037261517, 7.044143401826394, 6.679455991507198, 7.2237878747472095, 5.51102445900549, 5.888543938161388, 2.0046779660984217, 0.8477365476994126, 0.4641766646879885, 0.0]),
array([11.250052104485288, 10.318616591313864, 8.033969264724666, 7.680475843785107, 8.959724803327209, 7.865471176343402, 7.1269227075008015, 7.068475784917826, 7.226022149948345, 5.935124648086863, 3.9191365465829087, 4.610813377707057]),
array([10.714586922870895, 10.48922662023212, 10.630506661620128, 11.391963123163645]),
array([11.120239094318302, 0.0]),
array([8.106888639605177, 10.33038893514673, 10.804954922335819, 11.019993475328757, 10.032669010226806, 8.041905528469428, 9.994372116939862, 5.6174326726726065, 4.780447449130336, 4.591396506772413, 3.7725367976286446, 2.785629517896178, 2.1219571571007325, 1.0762882553047939, 0.9000012873175606, 0.9805476473071676, 0.5815074894736597, 0.0]),
array([9.188095050489771, 10.72684340042084, 8.472030695517983, 10.88485237531729, 11.06228467614498, 6.57741102488834, 6.016946942496981, 5.6307401447871825, 5.862702833120083, 6.915035591672252, 2.77508297614988, 3.342624966495271, 2.1953251775475446, 2.1853196198059166, 2.3430568659357007, 1.5250406666979768, 1.6350907990850694, 0.05376889491539469, 0.0]),
array([5.691721119071813, 2.771390798137091, 3.38424594345572, 2.3510775269318382, 1.5608291072173996, 0.0]),
array([10.494039289349729]),
array([8.878597736078001]),
array([10.216149341184252, 7.156495653631002, 6.420472224319657, 6.347153513428872, 3.8958018305817914, 5.146926495117768, 3.21702208260118, 3.5029615897257473, 2.7248052321071072, 1.4644249348446916, 0.7275694909380125, 0.1614934626971517, 0.0]),
array([6.994422576674813, 7.006522520773432, 6.8708796093348035]),
array([7.326482321564267, 9.807287385244804, 9.435289565197895, 9.281878527465954, 7.064460742301596, 7.072072215531385]),
array([6.296652793392615, 6.264096409799123, 5.604685957292067, 2.9328282407743966, 2.1636536166639635, 2.4445572459164957, 1.3569222738449775, 0.39865310950315197, 0.0]),
array([9.165110126061302, 8.003449785298365]),
array([2.94493065499636, 1.1316777153062336, 0.0]),
array([5.797516351616776, 5.037083111401285]),
array([8.507131162136004, 8.094330267320423, 9.191815984823194, 8.083177537577255, 8.721674529922883, 7.471942897365263, 8.304110012047623, 9.527648394229063, 7.032644531583705, 5.769650246977387, 6.466929646272444, 5.958760124427454, 6.622370484537292, 5.803551641179924, 6.0902936997973836, 5.753970890293981, 6.693812674417763, 3.8155519925046884, 3.4864217343326276, 2.783169829922025, 2.1036253419911595, 2.0749539396743217, 2.3989289083423455, 2.1050834218585273, 1.6755291635565521, 1.0402671037195153, 1.1685376835859291, 0.47301655339616017, 0.3575164273646488, 0.0]),
array([9.017329984192832, 9.015803802749593, 9.682821604701823, 7.5204695804637876, 5.7532626255052834, 6.485823419454089, 6.8426192131203925, 4.581210233469319, 4.552619534522435, 3.0646335572719616, 3.2630787385435687, 3.2081379151051004, 2.9383862628553543, 2.0217541298198825, 0.0]),
array([5.32857139065098, 0.0]),
array([9.049281227478547]),
array([9.11450808864583, 8.826115447288945, 7.284483827542463, 8.637385663507324]),
array([9.201196856398033, 8.289985959778367, 9.108801052233249, 5.761897713755298, 6.594910948585029, 6.979960170222389, 7.150448027417471, 6.230006255825371, 7.078049216502851, 5.742645754678751, 4.788588377013394, 3.192702731271345, 3.310991051247946, 3.561228196497441, 3.5780447324259153, 3.431412252145401, 1.8665421561485755, 1.8569919481499646, 2.502897755496379, 2.540150888550734, 2.0399090023992854, 2.18838047842986, 2.0459981382112744, 2.2396768524814905, 2.227324481529716, 1.4734013318823682, 0.16276759170880495, 0.0]),
array([8.600712858450427]),
array([5.7181254134204735, 5.6302872576849925, 5.898386717755, 6.728731739556474, 5.486527477852533]),
array([7.831642916837006, 5.538575854918873, 5.0411038859766135, 4.326231698116102, 3.4663934277594244, 0.0]),
array([5.594966345402653, 3.507759495786451]),
array([8.286616749651559, 7.134861291205846, 5.789077431302715, 6.747000649485887, 6.892127379526935, 6.0731831062057555]),
array([7.5693784913464555, 8.749354756063312, 7.350803905559993, 6.902304553547239, 6.544322772816079, 5.5750503457555105, 7.08603769391122, 4.149485812236123, 2.8108684726359416]),
array([8.92928577674922, 6.006043580583254, 1.6236074373822689]),
array([7.847912662430879, 5.4145359899153105, 4.555985678724441, 3.979609120219898, 2.6433744011033697, 2.0077643902285933, 0.7856952005584179, 1.652466534361369, 0.4362139010606488, 0.0]),
array([8.560202323557283, 7.332741497371181, 6.482259335908706, 5.459390167607362, 6.579037071942343, 6.459084484359005, 5.366104521517669, 4.366099661106201, 3.3389318446191556, 3.500527016959473, 3.1929743347161628, 2.9207785593122244, 2.630267282092937, 3.3659307666237663, 3.177694655052142, 1.958412865132409, 1.8581780140502016, 0.947716668219499, 0.9433362044951563, 0.9154911909501187, 1.5092971441402834, 0.31014885361284605, 0.0]),
array([8.519764700781531, 6.942010670242628]),
array([5.887650108028648, 4.826140751875047, 4.339451961800179, 3.6768598412177824, 2.61258140839705, 2.9705022443373954, 3.41427160671561, 2.033292807869401, 2.3770325821420744, 2.1735793135406443, 1.3263493984901993, 0.5618156172706907, 0.0]),
array([8.6980750551707, 8.325705648997227]),
array([6.068860307834237, 6.582863588269661, 5.761537951808334, 5.521513687850858, 5.755463181917251, 6.516441027519242, 7.226272429574598, 6.7662974387358545, 6.07089959322214, 4.487137177512833, 4.139953292452729, 5.175887380826875, 3.8606563078734526, 3.699984472978902, 2.9350855052211084, 2.8124282380859507, 3.580317078244093, 3.5940936983039724, 3.572679131185464, 3.31648372535497, 3.45837873003462, 2.5314797479464914, 2.0093852536766557, 1.9074431588844445, 2.1854783465752825, 2.312350804620176, 1.5114904884988367, 0.5873987508394449, 0.5939101282938765, 0.07921937692447117, 0.0]),
array([7.8857162414638236]),
array([7.9554877959766, 8.306066407452162, 8.236911288391502, 5.686965487092579, 6.781056047948719, 6.696649389320331, 6.898520111256718, 4.988272771622667]),
array([6.810087622789697, 6.8156746790512175, 0.0]),
array([6.307221607902827, 0.0]),
array([7.962972864380332, 5.513051189353787, 6.4435861676891015, 6.314585479964977, 6.992033911200954, 5.885711666952915, 7.066864037909706, 7.137407233006153, 4.079689974126347, 4.7066112815975805, 3.5349249235712623]),
array([5.43736303162704, 3.66426630198454, 3.3258142341742403, 3.175429334343197, 2.58577689072518, 2.094549137951074, 1.2142299148970823, 1.273018562705467, 1.247387766851949, 0.4952930055981195, 0.0]),
array([7.360372834355864]),
array([7.587055456124093, 7.567869642749768, 6.303397333315809, 4.386743202153822, 2.0916867332970077, 1.950079627510517, 2.3294111264744553, 0.06664695270234547, 0.0]),
array([7.096772851440793, 6.597441085134562, 6.909964313367793, 5.251037040115083, 3.1825631397890675, 2.9079933416821526, 2.5032756328355705, 1.8268801105313468, 1.5197219633701657, 0.08663747818661109, 0.022458740046816153, 0.0]),
array([7.874098822625932, 7.57065599764767, 6.847879283629746, 6.956674346364082, 6.775486463238868, 6.835620267442318, 6.9609541869404135, 6.743693625587893]),
array([6.051092547724262, 6.056783131980276, 5.751983406636608, 5.598370652461121, 6.774156123604596, 6.995757085033421, 3.338942694054052, 2.736295028098264, 2.9750632146530123, 2.622555510995678, 2.140175493295416, 1.9038724526816189, 2.052538757061924, 2.078905035159792, 2.0930496549559474, 0.7369502369196917, 0.20019191304739414, 0.7228342758399766, 0.12089164524206132, 0.0]),
array([7.324678024672157, 5.435108809937192, 5.84428669027772, 6.8877028257148485, 7.14020555793025, 4.460735136316952, 4.319098597356604]),
array([6.47274874254613, 6.1178275141687966, 6.826053504352763, 5.9701239171313105, 6.164310201608097, 6.965603210225904, 5.700764887667756]),
array([6.236209970804684, 5.964012341516711]),
array([5.599680296439346, 6.627975665233008, 6.462802604594139, 5.8106965694560495, 6.839842947520397, 5.390177221403629, 5.99695751697509, 5.7573995234275115, 6.233017076856445, 5.461567040724858, 5.458106304986346, 5.942711250875232, 4.649151308051317, 4.942853496452636, 4.722534323686533, 5.008167525878865, 3.745707979533346, 2.983751801559305, 2.8796309499815567, 2.709272278846634, 3.575935335431828, 2.8393813515357187, 3.1566821574628787, 3.439624433919785, 3.5929416546830595, 2.066780158933784, 2.2646555557458736, 2.0038885068550654, 1.9774153677869972, 1.9202741088148625, 2.0457748392638133, 2.0277177133867657, 2.114043477686243, 1.5530330616540946, 1.3721023084462582, 0.6485124184120303, 0.019879924296724222, 0.0]),
array([5.426403411994519, 1.9049848495098236, 0.20365116710231324, 0.0]),
array([6.7895274252962, 6.091208876192387, 5.459163055993552, 5.565152790212531, 5.37859822461101, 6.5965173299600695, 5.977579185554274, 5.715367279985233, 5.010260662917411, 4.039129600067216, 3.2522132039586658, 3.5140431041682243, 3.1934046601501853, 3.3574271712525925, 3.28079417753307, 3.1260248685366374, 3.0329067272687755, 2.187080787353829, 2.3410572814360595, 2.3479151683119506, 1.8845651811155522, 1.6583269017255224, 1.4098549064683064, 1.1532446940176624, 1.1230884619154557, 1.7973577307990938, 1.4913974618728325, 1.2865300338959116, 0.392059352932643, 0.3240797570209259, 0.6243780351573955, 0.0]),
array([6.287363071076848, 6.701626833445372, 6.695502024617899, 3.1388926080727875, 2.502313621831977, 1.4578308750514175, 0.0]),
array([7.073970737070292]),
array([5.8037120114561045, 6.145701497396945, 6.026290588564306, 6.663253064472114, 6.328431285508612, 2.615432181263225, 3.34984230148057, 1.8731056559896004, 0.0]),
array([6.704948136351251]),
array([6.2800056825698976, 6.467872407473231]),
array([5.5988343457017535, 2.7884205526741406, 2.9832893221850396, 2.6385829209152902, 1.0778542276077996, 0.7718285926040893, 0.3344427481769551, 0.05814823809485434, 0.0]),
array([6.316231145838197, 5.332026823656471, 6.03390953157886, 6.189921603612149, 5.632420088534366, 6.249225706707087, 5.359746511458164, 6.544792098768909, 5.870819013862039, 5.625237005723274, 6.81356518101726, 5.0257149968947, 3.9199921811436047, 3.870436957256132, 4.4124719447599965, 3.461461010226298, 2.6594610893931656, 3.574221858749721, 3.123724794873455, 2.6331580064837636, 2.661154213829793, 1.8596954195923088, 2.5118355065550135, 0.882901237633702, 1.0450477160525562, 1.2061373110524074, 1.7959922045981556, 0.8751835218809554, 0.2474279554918255, 0.40443561061725486, 0.27279735211630374, 0.7791563178393274, 0.11186968749386783, 0.0]),
array([6.14356774794644, 6.801219177627611, 4.257014969946138, 3.0651751347755107, 2.040300519048819, 1.5513439188459022, 0.2969730367053874, 0.432239657162166, 0.0]),
array([4.572496260413327, 2.89961400490904]),
array([6.471685468321866, 1.879099177006685, 1.2960683256554915, 0.9735046831807719, 0.7030564548773945, 0.0]),
array([5.813942081383603, 5.331839170907179, 6.0700280539170794, 3.7806346006056772, 4.347963657081139, 5.040072846602534, 3.4742342935206607, 2.5650670843430827, 1.673623876876224, 0.6873410445082886, 0.0]),
array([5.6421237272042655, 6.122579228509594, 5.770988530298114, 0.28311805910449495, 0.0]),
array([4.253134446646951, 3.4638561782595794, 3.1113492555762274, 2.4180037195863346, 1.624068652151845, 1.141868611103507, 0.9942239268161761, 0.0]),
array([6.001829450841677]),
array([5.716977787272275, 5.313823838322458, 4.478232972061323, 3.8673800771820748, 4.6232175002594325, 3.3559729234624225, 2.7453411111296955, 2.051977706673093, 1.8476139575805286, 2.442511501983099, 2.0757369827011978, 1.1326212453374271, 1.4936162078594324, 1.4400190542325686, 0.6193221624776883, 0.3213184551762831, 0.5058430074622247, 0.0]),
array([5.819472411424774, 5.939526468953681, 4.017358682103052, 4.930344959958264, 4.984517190533412, 2.6814209501743775, 3.511573088897505, 3.333968975121017, 2.2928855930017336, 2.1349475644016462, 2.537704721870058, 1.2825872952427821, 0.6504203980287901, 0.0]),
array([5.553630508319499, 5.646793719425468, 4.756956991715728, 3.762178130765238]),
array([2.752727215125493, 0.9032793632988111, 0.0]),
array([3.0616112178341814, 2.5122557312272096, 0.0]),
array([5.614424861436868, 5.387916333690675]),
array([3.527070923542372, 3.0835307325741157, 2.251130544901657, 1.081668864504706, 0.8684178127523855, 1.0174481741127614, 0.0]),
array([4.288544019443455, 3.286767347791893]),
array([5.340232254329976, 2.9006158169433514, 3.301078558653826, 3.503377735468701, 3.588336072983157, 3.0405990784490626, 2.8089833257850407, 2.5226629221295656, 1.9546558878561573, 2.2852967582116612, 1.5100612343480546, 0.8945758356184782, 1.0134132139927157, 0.9694161984012895, 1.2356129175411659, 1.767219413433621, 0.5439160307134633, 0.48639734550332037, 0.0]),
array([2.2037101478324996, 1.1261332769894994, 0.0]),
array([2.613600987486233, 3.442428109087122, 3.4853707658350315, 1.964745995434054, 2.450661906206425, 0.7870329440669268, 1.052100001900903, 1.3562986087661124, 0.956919827182838, 1.2848028474936868, 1.401962992061584, 1.6924889692271763, 0.5570394534777112, 0.6707886453268596, 0.23684158156513124, 0.0]),
array([3.015816278008153, 0.0]),
array([4.380078256245629, 4.967884962138173, 4.929846980978166, 3.562192716370473, 2.873355120535277, 2.6589013183705488, 3.072331667928772, 2.5726410770904073, 2.5494776663250236, 1.2716008023588048, 1.343403664584641, 1.3860435346832274, 0.0]),
array([3.650745193451921, 3.9807804959739475]),
array([3.314829232818319, 3.0655803470086713, 2.2471322412831323, 2.45360385742157, 1.5547896964561188, 1.0696461534110564, 0.0]),
array([4.092416175824296, 2.6296512601606667, 2.978842459816918, 2.404971536272935, 2.4097374509495744, 1.0569419602944876, 1.3715289620931286, 1.2956519191485942, 0.0]),
array([3.1404994079144517, 0.0]),
array([0.8079904838195374, 0.0]),
array([2.0192337151607695, 1.9011733765857253, 2.0064082814620816, 2.475714797445653, 1.5436868426955883, 0.0]),
array([3.6321077722508326, 3.063004599093257, 3.167318828059776, 3.0884355179565013, 3.5424087806744886, 2.1773879828445057, 1.6262052724345084, 0.8009398692396228, 1.0488007316401116, 0.9278127641923898, 0.04494335729670505, 0.08811265053233189, 0.0]),
array([3.78838234318044, 3.888696696211681, 2.869985608294879, 2.997529914508178, 3.102024853753817, 3.3415404147073295, 2.3972176078649086, 1.1312017894951663, 0.9109651065647222, 1.2882247254739043, 1.2993985953602087, 1.2270528339042321, 1.559729022099164, 0.3595739673626474, 0.09096166511371778, 0.0]),
array([4.190575341605461, 3.4468632227498865, 3.1109246667147925, 0.0]),
array([3.8211646141050757, 3.010271651374787, 3.0384456319391657, 1.870769252741524, 2.2148013743096175, 1.6423967240180137, 0.0]),
array([4.320243354984773, 2.0274788761188445, 1.8583944868265792, 2.5115119906501877, 1.1051355385503903, 1.3777721627787574, 0.3428438337524634, 0.0]),
array([2.938986367032842, 1.8570400410138204, 2.4154269090991165, 2.154731590229792, 0.6346292009964734, 0.0]),
array([4.1242972881251845, 3.210439773180993, 3.4429484540436373, 2.9306710707393164, 2.7049422173350957, 2.0385466651278525, 1.6368758869713556, 1.1537641809550454, 1.5074996409372492, 0.6021177178705326, 0.0]),
array([4.348962557530599, 4.508397868390042, 3.1613370394288838, 2.9839200438354525, 2.2823549078750145, 2.3405788619604864, 2.1616738659727313, 2.0168817617317405, 1.1411432411479918, 0.0]),
array([3.923791852359153, 3.907065498709781, 3.46299987501661, 3.0609728566342493, 3.573875406489889, 3.436025672588627, 2.9468764029701364, 2.258400886144057, 2.363648351890902, 2.3075034791379867, 0.0]),
array([3.3983325111162457, 1.5328113885657755, 0.0]),
array([4.109288488365012, 3.4710202970270188, 2.615758038545104, 2.8638526385341914, 2.923801487619347, 2.1736255326661946, 2.461299091211878]),
array([2.6607544974028476, 3.5776236173794174, 3.3556732136110505, 3.465155282427887, 2.948846030522774, 2.5767149323855474, 2.450916522119205, 2.353643445582325, 0.8337285257085902, 1.2048649637835394, 0.6830170838484477, 0.013222508374246975, 0.0]),
array([2.6446696980805062, 2.7463699481222736, 2.376493054223541, 1.3765055142051459, 0.03055975913939374, 0.0]),
array([4.050722782315832, 2.936096764982744, 2.6425178721484097, 1.9960531277625595, 2.404179301343619, 1.9618500912459116, 1.4623834750458589, 1.3516354149806045, 0.0]),
array([3.0605781466786417, 3.187464250270903, 2.9210620060670975, 3.355510944186794, 3.3383195564012698, 2.8355295684446706, 1.8408969853775168, 2.0246034708846095, 2.3055432122781125, 2.4726745039366236, 2.0090397447108925, 1.762861726428296, 0.6848332836568772, 0.6818810473462804, 0.0]),
array([2.872126773721258, 3.3799471662584675, 3.2654026121890967, 2.5303935036616054, 1.8574466804953431, 1.2106439406547302, 1.322085577287592, 1.386569894463579, 1.6072667750288685, 0.98669000830511, 0.5585632078675509, 0.0]),
array([3.7502429495713456, 3.350449185429003, 2.2449210407738835, 2.246779963945827, 2.1596216442168936, 0.0]),
array([3.920191434977523, 3.6135659009023833, 3.429716247196005, 2.937760665783371, 3.4425391950694655, 2.6370603155068224, 2.2889564249381693, 2.5724598693669676, 2.545694612602255, 2.0277515956081844, 1.2714193946439694, 0.0]),
array([3.3891670810256587, 3.1354256619485144, 0.0]),
array([3.13160480807951, 3.0259073374098984, 3.3489516751730823, 2.960925910371382, 2.044603837743901, 1.9223010118129893, 2.0756165121924006, 1.9379337096615807, 0.8478164249745224, 0.3380230754195889, 0.5757591297998285, 0.014313897946597584, 0.0]),
array([0.16910407336077693, 0.0]),
array([2.8988017871119216, 3.0916554013358284, 3.4699787341806356, 2.8717283893885934, 2.901614401018041, 3.263936434308211, 3.4646191309310717, 3.211735849579917, 2.601342446881486, 3.575023114090156, 2.95269606123274, 3.476403776280342, 2.95355934136621, 3.1358288416015077, 3.5007282622646003, 2.948993704979714, 3.433979389132148, 3.5783398514879505, 2.843746562811471, 2.397951592742158, 2.323975355534708, 2.328680329177421, 2.045368851736804, 2.495223412150051, 2.0159229939752343, 2.547785731904544, 2.348763432319131, 0.8402032491154945, 1.0584550936469856, 1.567479520751886, 1.4943148446513472, 0.6451118402217979, 0.5736017271411458, 0.22285375109410788, 0.16627846485988473, 0.0]),
array([2.852665470895769, 1.9794305038898994, 1.5720473912387722, 0.0]),
array([2.776366830491144, 0.7072758675881415, 0.0]),
array([2.901542817883535, 2.886631930849572, 3.0437587658057024]),
array([3.399859008470766, 2.919649554098956, 2.145608455947808, 1.9594057507866613, 1.8983236893889377, 2.1207044486775315, 1.994282948989703, 0.7347934792286868, 0.1368086300387069, 0.0]),
array([2.3595443814401755, 2.500490772241221, 2.299530722814771, 2.47101637913524, 1.91963489832822, 2.4087042359127957, 0.2264541184908052, 0.0]),
array([3.1815325398673697]),
array([2.8411959930981174, 3.280607710395681, 2.4007760739601434, 2.2427064247336608, 2.1309214989364453, 1.1622069692794605, 1.6666761496332017, 1.5755414404718564, 0.10177120254551991, 0.024567749329035965, 0.0]),
array([3.004042685658049, 2.8129485245693138, 2.79230832058935, 1.6501295915486918, 0.7963989824568365, 1.7535335922892186, 0.0]),
array([2.708538770460786, 2.451775935722826, 1.6328145813091017, 0.0]),
array([2.756962523988593, 2.967699121474844]),
array([2.4055329411735418, 1.460248776475308, 0.3224113591121313]),
array([2.973034221304373, 2.9717768931666324, 1.1854687117860077, 0.9822878160384229, 0.0]),
array([0.39870611390333377, 0.0]),
array([0.12070241625165434, 0.0]),
array([1.3717821466660203, 0.9839854437324396, 1.6401310086174807, 0.0]),
array([1.9722512705785, 0.0]),
array([1.1978358698391838, 1.1099723017659047, 1.4746420084335479, 1.056761301379129, 1.1464073310410514, 0.6943295739772631, 0.3210434571942749, 0.7470189155891086, 0.0]),
array([1.8462358251817785, 1.7659736881942143, 1.2377605988679985, 1.3045166603385894, 0.2874778015549801, 0.047207055146081936, 0.0]),
array([1.9418756894266234, 0.0]),
array([0.9447782250685799, 0.0]),
array([1.9226226036582468]),
array([1.417685659230962, 1.2002778626694084, 0.7167526405759107, 0.44169616643250276, 0.30140422218966495, 0.07846474235203921, 0.0]),
array([1.0727642444246925, 0.03443012877897586, 0.0]),
array([0.9533589544744132, 1.2705597540532028, 0.8667356569563992, 0.6107109877271292, 0.0]),
array([0.8545850631829425, 0.2320356482850945, 0.0]),
array([0.5900224336849373, 0.0]),
array([0.11569097359333518, 0.0]),
array([0.9043258386033296, 0.0]),
array([0.34697775790370056, 0.07612319794377975, 0.0]),
array([0.89135296750663, 0.4407561718507048]),
array([0.8494021467180866, 0.837076430514188, 0.42657277783883607, 0.05109317180072524, 0.0]),
array([0.9099005428743564, 0.6776706085534528, 0.6843150504952462, 0.0]),
array([0.592089312512991, 0.6239468737776394, 0.728846088583132, 0.09651446924074814, 0.0]),
array([0.3829078824514771, 0.04180060270434223, 0.0]),
array([0.7941767106947362, 0.12254679473309205, 0.0]),
array([0.31811086504095243, 0.0])
]
d = [data_1]
names = ["59"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T13', 'T14', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T32', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T48', 'T49', 'T50', 'T51', 'T52', 'T55', 'T56', 'T57', 'T58', 'T59', 'T61', 'T63', 'T64', 'T65', 'T66', 'T68', 'T69', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T84', 'T85', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T94', 'T95', 'T97', 'T98', 'T99', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T109', 'T112', 'T114', 'T115', 'T116', 'T117', 'T118', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T127', 'T128', 'T131', 'T132', 'T134', 'T135', 'T136', 'T137', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T153', 'T154', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'T165', 'T167', 'T168', 'T169', 'T171', 'T172', 'T174', 'T175', 'T176', 'T177', 'T178', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T195', 'T196', 'T197', 'T199', 'T202', 'T204', 'T205', 'T206', 'T209', 'T211', 'T212', 'T215', 'T216', 'T217', 'T218', 'T219', 'T221', 'T223', 'T224', 'T226', 'T228', 'T230', 'T233', 'T234', 'T235', 'T238', 'T239', 'T240', 'T242']
def get_taxa_names(): return taxa_names