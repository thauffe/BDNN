#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.81387051367059, 31.098268407313142, 32.039047796521075, 29.026519939253472, 32.650783743891964, 33.61473539012395, 30.482592316511557, 29.109090373037944, 31.092948391899846, 32.38960696085956, 31.001062872777574, 31.544537719934205, 29.7233638553655, 32.896425036319, 34.0232342255874, 31.635493767760615, 27.182238493263903, 27.90930575786108]),
array([33.11690226426852, 31.072671322888496, 32.65144397972284, 31.83315664858343, 31.881517907999623, 30.49883980301057, 31.48041915124993, 29.065102329994403, 31.662897544250065, 29.39074906822945, 30.885194524558184, 29.938618847702287, 31.982792863173934, 30.636364060033383, 28.387374234965858, 32.93201619484172, 29.840787431256146, 32.80233430855461, 29.178990193565017, 29.580892926070447, 30.486753946100013, 32.60862926368869, 30.037538291306845, 28.132244700941484, 29.152931249826366, 32.18089707076876, 32.234203370702005, 28.754044763674457, 30.106240527317386, 29.652513146444655, 32.49764727939048, 29.184978085964964, 29.83458680661575, 32.23663911171329, 28.122360892517072, 29.797452931741077, 32.179614148820676, 28.31147565859391, 30.282066561122708, 31.91914417597935, 31.149775064093568, 29.829546916512648, 31.349384540294977, 29.43469374356852, 32.833081206155526, 29.62533913350373, 29.82578496241114, 32.3401493771064, 29.80603257434056, 31.90028678323082, 32.4093457648478, 30.947063017200705, 32.79554968752775, 32.15585285129828, 29.76156966971622, 32.86154126307858, 31.914291991181464, 29.38694131072981, 31.67644310220569, 29.87644354447192, 30.9892905017673, 30.68838665228056, 31.280835414979116, 30.21773215059082, 33.055836065100074, 32.18367201950329, 31.879402288431347, 32.413062859239666, 30.36289752495685, 28.6151167475534, 31.13023147401149, 28.36141327186406, 29.152363786455012, 27.93394786890116, 23.699089109736256, 25.76932107818114, 24.445044734271264, 27.658986322391694, 23.582685367368807, 25.78565773082305, 26.232990405043292, 25.94543738468576, 24.289321285014495, 25.25124591610343, 24.60177909327622, 26.00036178317515, 28.01349139513903, 26.16214519864632, 27.80842430027902, 23.986508963821404, 27.2399079378422, 26.417220510359083, 24.78117031228024, 24.06648488974379, 24.440132955561797, 26.362803451445515, 24.225573338307377, 25.980362348170132, 26.88091393395357, 24.10612121879043, 24.607326712454665, 27.178773446463424, 22.319928236649215]),
array([30.424447363954584, 29.853293245525897, 28.756997784273036, 28.373999573192865, 28.583294249934966, 31.278405291607502, 28.657316712485144, 28.980615902312188, 28.38474216435258, 28.99893292099839, 30.87531504908076, 30.186061440019685, 29.106092722356415, 28.356288768175123, 25.248689802913464, 26.733218404919107, 27.079245714997324, 27.71044335633244, 27.9910791494648, 23.038070277286753, 21.003168988890707, 20.99012585113808, 20.801993428722056, 17.884574733119507, 17.467124219842102, 15.393280561668421, 11.765670914076017, 12.123671531677319, 13.697420417770427]),
array([29.874914982083844, 28.71855030221667, 29.60060193695402, 28.85108239415056, 30.263166178234798]),
array([28.499601890499367, 29.2109051509494, 28.324232055166508, 28.703046655891523, 28.39315776618918, 28.5591605963813, 28.835213640730817, 28.501833121817988, 28.513253778671324, 29.139787232104045, 28.569180027174507, 28.107095249940123, 26.904103012086942, 25.675789217859865, 23.035095748641744, 23.603316673523942, 25.389420763710614, 26.45242430695676, 25.121715805869545, 26.29369525097308, 25.930554648264955, 27.84714726337823, 24.017435853634286, 27.602223146101185, 25.631888447208805, 22.66491266224147, 22.952013216112746]),
array([28.62885371599567, 28.25409997398479, 28.199218586206722, 29.062559435203788, 28.76825569113611, 26.410935117077294, 25.816253078377084, 28.00918823945803]),
array([28.27876308219622, 28.48940014834618, 27.796419062037394, 23.521999892716003, 24.399636514809853, 26.276366003089553, 25.135422141822765, 26.260671470572422, 27.398442606030862, 27.340833220066116, 27.642177367313238, 24.348573103496214, 27.612200367486825, 26.980690704636306, 26.348117722352065, 21.22509763828755, 17.32308643826193, 19.208330331990226, 14.152454721673926, 15.827806432288984, 14.761548575972842, 15.126862364261575, 12.479663842067694, 13.060316398515088, 12.330449021801037, 13.250874022666451, 11.805946187090958, 12.004171031414979]),
array([28.551900617891093, 28.565652451837753, 26.872151784802966, 25.668589226335552, 28.000758534899756, 25.803922451769957, 25.552323558140326, 23.845056086293724, 27.74159545437025]),
array([28.199011205268548, 28.1104707751371, 27.50457514287908, 26.0795219665203, 25.515938443300016, 27.541382684759693, 23.56024165740086, 26.533787770797385, 23.869658470316864, 25.010268489766936, 23.663417227623874, 26.05793621209831, 25.07719261114981, 24.670304568801782, 23.20352382279616, 24.33322253046252, 27.85945114111161, 24.23192643059567, 25.3388340405761, 27.169860701507286, 21.961171577444546, 21.649596796450837, 21.454057682596424]),
array([23.619213877812054, 23.111038373493393, 24.17139579296304, 24.533417173133724, 20.800999120982183, 22.72943873524086, 21.761868740853043]),
array([25.368496842112382, 24.813179284281205, 23.209011128422354, 25.4209811469814, 24.716254074564187, 24.52378570859814, 22.968034722771833, 22.8629220209286, 22.558691993117964]),
array([23.207164874199513, 24.2253715654764, 21.928819532715774, 22.543943035407718, 19.399973521878717, 20.035310488399666, 19.851247479647117, 19.156062507849867]),
array([24.23155905363553, 22.10178966367824, 20.66668980940021, 17.47009542232678, 18.812214126829364]),
array([24.586206921212625, 23.244034364030526, 25.130449821151515, 22.107033532252853, 22.73484705338138, 21.412933598236947, 19.631086850245083]),
array([25.34757877879234, 24.097158347837656, 24.74000418451726, 24.282369103090627, 24.977698897531877, 25.271336948574493, 21.302028470643652, 17.704651003824495, 19.059026904954276]),
array([23.891114559803025, 24.22156371005486, 22.310549216631838, 22.516498137324742, 21.027406174185757, 20.451262181875936, 14.64613170787133, 15.072934122121717, 13.119928078480342]),
array([23.337873713752465, 20.774582745325258, 22.93025062445541, 18.771139364223764, 18.827957807401255, 15.760492684248558, 15.78682448054679, 15.549395033956717, 15.640071874433463]),
array([24.081351716483127, 23.7721082904928]),
array([23.455627510070244, 23.7215183567449, 23.107616454935886, 23.914764371706703, 23.498332572121548, 21.735254552915098, 20.785179795925924, 22.76076213408521, 21.801904102315277, 22.421408859331464, 21.65810243505763, 22.193445833464118, 21.153394338108573, 22.863712851171595, 21.325707687269905, 20.419492404625828, 18.396680343451077, 19.07355903805699, 18.179707815262503]),
array([20.960652420059166]),
array([22.154468190556088, 22.458896990223277]),
array([23.174296430035206, 22.658758291059012, 17.794369392817217, 16.397564276149236]),
array([22.141453910218374, 21.00526315664462, 22.693651697291767, 15.789724711326336, 14.545696719252257, 15.684365249646046, 14.856381445184768, 13.68222166407033]),
array([23.5210842429793, 23.462582980109964, 21.932995885392867, 22.588625112470403, 16.49569851678017, 15.793153254113298]),
array([22.84189942922744, 22.41659673784988, 13.31599153441951, 13.318535156285629]),
array([19.013187081964816, 18.670941433498, 16.36382273669544]),
array([23.393302031146145, 23.23534783297245, 23.223908010736153, 23.47998277921116, 21.210019468465113, 22.683751147937556, 20.546162151081262, 21.08086521962181, 22.670705967798487, 18.720733257580576, 18.60242144009556, 18.740907095063264, 16.1508253177661, 16.175693034362148, 14.308051836699173, 15.093940197213605, 14.16390178244329, 14.613235283735253, 14.117728337322138, 14.89939859936556, 13.177592821537491, 12.225953429261281, 13.110034351126608, 12.485501388556093, 13.562753698111546, 11.821350952100255, 12.572820146457316, 13.539797625807257, 11.257704685412852]),
array([21.90159629317446, 20.33838266487899]),
array([21.269045618892854, 22.8837404784696, 22.403407303565384, 21.906010135958674, 21.950818881036582, 18.645807981612567, 20.061074082368084, 19.399429901930567, 16.502036222101914, 16.702992092399203, 17.71755849243644, 15.20110661641459, 15.559923244009694, 15.076601718173658, 14.963440262581312, 14.024607748221593, 14.93243470530972]),
array([22.489823844045617, 22.17783599758618, 19.526278179754303, 18.626278774442135, 17.23712735340586, 16.348552663580758, 18.318287240832785, 14.394111693832196, 14.024567284368992, 14.203497875844722, 15.505740194390393, 15.122973457231065, 11.923306123955102, 12.4582492393558, 13.6436216284005, 12.263752513997431, 10.64816226048321, 11.224816162118692]),
array([22.991787688910016, 21.738327040998104]),
array([23.062093105575894, 23.05471920751405, 22.860125556181874, 20.55947051893989, 21.85364994435779, 22.285626201640813, 21.54610540495055, 20.60355100506488, 21.342412721623408, 20.82076133719369, 19.124806243305212]),
array([22.268412519626793, 16.106862927138895, 17.634643360601554, 16.914356702074972, 17.550602972969088]),
array([20.69433107521633, 22.114621167758823, 21.77865615574427, 21.979184914099257, 20.836117589273837, 20.49945372743553, 21.77936573974192, 20.584682733960218, 16.955272810228823, 19.291788393462458, 14.748297603306515, 13.874468026020065, 14.632261912256977, 12.659581745975453, 13.225754517359846, 10.622732787514858]),
array([15.966930148211219]),
array([22.375651121917276, 21.616913273086656, 21.073369060963415, 20.525341632693802]),
array([22.448400431130356]),
array([16.785071925201134, 17.935610270616277, 15.108577764756985, 14.94409249071511, 15.080963020274325, 13.72057188094485, 12.433088327640206, 7.299153368308052, 4.582219692004433, 3.8884904100173525, 5.21091863923148, 3.0486241224387105, 2.4349888639198545, 2.230434077214725, 2.5009502421905876]),
array([17.75453976573501, 20.004292458100416, 14.995681169235482, 15.03968203747734, 15.351338294733836, 15.452157979855864]),
array([20.81927545317809, 21.77630743542603, 22.182975949999797, 22.0580071899431, 19.4909400065172, 14.366167537727605, 14.641075629295479, 13.925821141411157, 14.347977122549565]),
array([17.947230059548453, 18.597857096089847, 18.890732234080392, 15.139497803804959]),
array([21.034174009782895, 20.816169343370287, 17.67678689967212, 17.963447711867545, 19.12743200186237, 14.238775421981483, 13.13962832007473, 11.868368659535848, 13.525187280130455]),
array([21.458142208495218, 20.80656857874553, 21.72309261256718, 18.76822549057762, 16.574463400966646, 16.948663493099826, 19.270548399311423, 17.6248977999337, 20.06253562625667, 15.788003654739072, 15.319744243216494, 14.185950387322617, 15.403531660525998, 14.587193726036121, 14.05089281317328, 14.950767793732645, 15.217649884170672, 14.582596401047544, 15.759928485540156, 14.938079244940305, 13.465129045563453, 13.317605611993136, 11.654497392557701, 12.187314597277041, 11.970336135178417, 12.050046007105761, 12.968834898537352, 12.004591109347636, 13.817356461576688, 9.74771258457298, 10.598502845890467, 10.556844094846777, 11.246469014972952, 10.349158039752808]),
array([21.98915640843146]),
array([21.612160493682868, 18.818033084218182, 16.377268910878673, 14.95046778888272, 15.185099155826457, 15.454817165010304]),
array([21.014868308750145, 21.10938939809926]),
array([21.19859260732844, 21.345444662241878, 20.853233022317028, 18.298327611455708, 20.020158407783946, 19.765921409734865, 16.66220761622419, 15.047152016674147, 14.98644798042714, 14.539299246576363, 11.84531376482723, 11.773393716486659, 12.718507754619289, 12.152249160559625, 13.626101661155833, 11.03013366846437]),
array([21.025259263622157, 21.88289518800835, 18.31368075347235, 20.3867473324019, 19.236778308461446, 19.208088870744124, 19.128260536970398]),
array([20.99362036168033]),
array([20.860438090860786, 20.281814516496524, 19.118895822937773, 18.454870901186386, 20.15984341185493, 14.37052847650014, 15.802618490159706, 14.222736826869031, 15.160174519720238, 12.281005126025018, 12.402289977974712, 11.815874617206692, 9.201271338917199, 10.438878058885095, 9.495437246665603]),
array([20.768404556830433, 21.15182792240282, 15.394630543096087, 13.831889572747375]),
array([21.375706051854397, 21.79738477868066, 18.291133760636097, 16.88491729714742]),
array([21.6101064985459, 14.286015258080797, 15.741527477641736]),
array([21.224737946605476, 19.848909665694418, 14.607133741116236, 15.443359716516605, 14.871931928151378, 13.242607404705623]),
array([17.756257858768183, 17.11027467413328, 16.83165441543138, 18.61620046665802, 15.056801818514886, 14.84032134527808]),
array([21.435480215463723]),
array([20.099326843207844, 17.688353062650396]),
array([20.76936124552997, 19.60842760859342, 17.925246096379183]),
array([20.580225410871332, 18.967567408051817]),
array([19.15616478021621, 17.854712201738785, 15.73628393793011, 15.347238281674102, 12.932546151127001, 12.423600557468147]),
array([19.891074721025962, 20.329500222885727, 16.992463762461643, 17.158751175519754, 15.663558525459289, 15.06631587500855, 15.31185583619375, 15.22206441407598, 13.317465757590709]),
array([19.82171026150306, 16.636689911799053, 20.021532945767945, 19.617549656882076, 15.429773443478542, 14.059157840088021, 15.805220135315786, 14.0331297174483]),
array([15.691386040659696]),
array([18.151303410791893, 19.488440050519213, 12.294635722351602, 12.307817823965962, 12.35861210983684]),
array([20.719949698180674, 20.159822970891547, 17.690387138837227, 17.667406176875936]),
array([16.92545485625097]),
array([17.395175078920662, 19.66246263584697, 17.94725969918498, 15.298882901698809, 15.030511509102709, 15.543738470635006, 14.578473813459116, 14.886684253105734]),
array([20.09876961185985, 13.980061326886956, 15.163552797578728, 14.894713659907035, 15.381906809638187, 14.796247186324313, 14.609090826272967, 14.480049591431413, 12.812053561245534, 11.78479260987492, 13.141065015330046, 12.818015886055708, 12.931869475424934, 13.65402186291379]),
array([16.84163194655653, 15.560060886431662, 15.708342527119171, 15.76677557504405]),
array([18.098711206655278, 16.189670446346163, 18.73416938582219, 15.713746798516821, 12.381030693739831, 12.85513673839292, 13.57048477583163, 13.378593726588795, 11.606398114959607, 11.364863619959939]),
array([19.848068514497086, 17.34451429668448, 19.74041208184585, 19.161401739385983, 15.475328114187064, 14.492864441761796]),
array([18.349067241754373, 20.34610410945898, 16.086406488138216, 18.61224813043837, 19.245126963569056, 14.211637107653825, 15.487599974863556, 14.758585842925044, 15.107036214524781, 14.816431301548624, 15.844492998455518, 11.786326704004086, 12.97686042198769, 12.760692984760597, 11.910552113532168, 12.128437860862256, 11.98137409491369, 13.75634684509132, 13.180306464580308]),
array([18.28434743367369, 18.785230284679233, 19.3513323420279, 17.9342768197667, 18.86910343214928, 14.805177477115262, 15.960365919832771, 14.473471578287475, 15.708718940678887, 15.231224515053704, 15.356077260938216, 15.165338512026462, 14.755290385921283, 12.507025276763585, 12.867935705959077, 13.000331586625801, 11.681720579579181, 13.69408953162693, 12.252984968414669, 12.725213323184008, 12.490281214010897, 12.17405553033579, 13.776691544959691, 13.339265501275884, 11.46096321229876]),
array([19.158635672706744, 15.556745730398264, 14.767073669324132, 12.799962927119426, 7.807488639868014, 9.025420058576303, 10.800799402676219, 10.123950515540718, 6.311165457353866, 5.9729563208232666]),
array([16.54061983992353, 16.76948095615676, 18.739118562772667, 15.872725993963346, 15.706468659123015, 15.559054659553844, 13.64683802234383, 11.837135192476753, 12.426391557981049, 12.558428039061631]),
array([17.85332280377237, 16.49384010898129, 16.617077028846627, 16.46002034919373, 17.47346087853932, 14.555410997729258]),
array([17.88468584265302, 18.347013558725862, 16.399700516386343, 16.0200989202609, 14.597658526292259, 14.78515955956873, 13.896688041251895, 14.081745659601797, 13.656351868503327]),
array([19.522254196357576, 19.037309603255835, 17.803028110551196, 16.902320116064214, 16.83641988476273, 19.623269236103187, 17.519700249889258, 14.942728937621906, 15.716258663600176, 13.861144193494612, 14.161039296153001, 15.843252128824222, 13.918376735253265, 15.829711657022408, 13.940322856003107, 13.443320230937829, 11.831470726334308, 13.508482030303522, 12.092588635609326, 12.269293040169616, 13.431907298227099, 13.701492915219736, 12.035847812409452, 12.60732972372242, 13.405358405031393]),
array([16.937723848343147, 19.596870186052623, 15.21037883545692, 14.894218767958183, 14.363540501711018, 15.194135965423245, 13.87161318597885, 13.559445190048233]),
array([15.919148622133454]),
array([16.548253954826304, 12.536667168283364]),
array([16.632183254209245, 17.4383240468605, 16.61918810826311, 18.34110450799497, 15.039796123479396, 15.905308271683623, 15.701172084139184, 12.176177390495281, 12.38780558955054, 12.118913379524736, 13.499020711796213, 13.36207450620889]),
array([17.188615233170783, 18.130065117297747, 16.711782737746468, 16.86075560178272, 16.916340180021443, 14.96872404713283, 15.91764334276805, 15.517092355075343, 15.142698492421594, 11.823420079051512, 12.730290571667537, 12.920460823436933]),
array([18.511068141780832, 19.203498654570815, 18.837669093201303, 16.232085871797423, 17.3922491972123, 17.449192515993847, 18.076258825648583, 15.67238011789569, 14.552856476680263, 14.42420017982424]),
array([17.371854281854926, 18.68542330848031, 19.56337373932917, 19.444916306895173]),
array([17.146864335444658, 17.182353113621318, 16.332698853341917, 17.577060377136956, 17.738016153421942, 15.454649345064993, 13.850724422497183, 14.85603267318936, 14.137417003273004, 14.227326882902696, 15.733833833437682, 14.621276329674345, 15.95487780367149, 13.69595141243925]),
array([17.47270304843775, 19.425281243157063, 18.507568773013396, 16.953049526275702, 17.994406917326902, 17.512619109917004, 18.70846929957012, 19.099810749520167, 14.20075670885355, 15.19523306343687, 14.475304886397506, 15.876355084479124, 15.665019524614449, 14.1400953585795, 15.04259054739983, 14.069457640129949, 15.246705426060943, 14.059430897387196, 13.310677879674921, 13.402302475186584]),
array([18.317288185248223, 17.069845710540406, 15.317409840833147, 15.324355664048724, 13.27060075996474]),
array([19.26913605451127, 15.934080484196592, 15.93249473009791]),
array([17.295651534805042, 16.49873154196444, 16.242870528519305, 15.860505270477407, 15.943268760007712]),
array([18.604407506760936, 11.897853368043938]),
array([17.47742690188484, 14.87758625701328, 14.125118590305402, 15.798783758933537, 12.839433100460717]),
array([19.11183295535505, 17.100276954358897]),
array([16.538131979465028, 18.311549052059114]),
array([18.246838166567453, 16.70572797487647, 18.11141142051508]),
array([15.48617516861685, 15.88276063259864, 15.125578652267944, 13.27696812783375]),
array([18.470250031076358, 18.94428770480018, 18.534633559074848]),
array([16.1767882315124, 18.638083987397962, 16.598764713373924, 15.278304060705773, 15.005334960146248, 15.871402578020511, 14.048213405629047, 12.185126336410292, 13.032955670826068, 12.939474928073881, 12.770356881450574, 12.671320603785766, 11.597236846700106, 9.986124707568376]),
array([18.433737110028698, 15.566153354315722, 14.322144915714619, 13.584160049748824]),
array([15.54394872594611, 15.617186183383756]),
array([16.848720875202968]),
array([18.49073363426493, 16.95427736498127, 18.380298576538028, 15.271982007545464, 14.223742847576265, 14.66330195376042, 13.99578333084621, 15.764283674621725, 15.100009647451245, 13.179158204096424, 13.511164157464082]),
array([16.60904448163464, 18.16433012054408, 17.04960294254024, 16.97248014040835, 15.926504035485674, 15.439003466644511]),
array([16.207690311012563, 17.03506105360798, 18.34789514436139, 17.616332633076652, 16.811846406393677, 14.783991680691049, 12.922680274177466, 11.853104972585598, 13.495706793731046, 13.594343651674633, 12.306878460779988, 12.18157281589955, 12.38638420532815, 13.663552069676332, 12.658696827226116, 13.685132472553905, 13.42202230297122, 12.055800455042897, 13.409192567776643, 11.621810054567899, 11.50499533758162]),
array([16.861951817748743, 14.66818063983686]),
array([15.89311354829977]),
array([15.068239549974859, 14.112030834614254, 15.73672982264889, 14.337309142776766, 15.408998700732148, 15.118693138187764, 15.827740118225542, 13.088566523633641, 11.693675526785931, 13.712014554725144, 12.605126222839433, 11.749164863956738, 11.954961715638744]),
array([14.960436225100564]),
array([14.45950978429963, 13.95162860610285, 15.42343202953592]),
array([17.58311257946996, 17.274504757313526, 15.664219896182018, 12.800812985469298, 12.437758936507908, 12.888581225046693, 13.183432559902514, 8.608575579762293]),
array([17.632456731939236]),
array([16.65182865784477, 16.98643300394248, 15.278675186224667, 15.248595070113094, 11.973377785128774, 13.19946801094949, 12.598504448130784, 12.229909926658097, 10.976630926329968, 7.877204899334146, 8.547316700646748, 6.3145384513901845, 5.78350048483923, 5.353516323523203, 6.420282435012526, 2.5964460286294653, 3.575933990212933, 3.3418482331371635, 2.283425312493476, 2.548220583770721, 2.3425811979128817, 1.9982723485486367, 2.3991851218277924, 1.0011603798116686, 0.9182453230876776, 0.5359869354880523, 0.43853965889451857, 0.13099832605029493, 0.6288988910869429, 0.34871456610242446, 0.22476784960064955, 0.0]),
array([16.40177396917543, 15.35962746648849, 13.869080427797131, 14.801990579004276, 15.519596995987945, 12.052716049175931, 12.323399610490725, 11.957802943219821, 13.383543580945382, 13.465667193889582]),
array([15.823669729019343, 15.38055184685647, 13.07068845309033, 13.104098111688101, 13.223884113029694, 9.40456487535684, 8.061957556741746, 8.729185947783533, 8.176201076819492]),
array([15.691757574534167, 14.913300205669719, 12.302699150716593, 11.526578553842196]),
array([14.628587430754957, 8.940466048705513]),
array([14.090831830128483, 13.310217192130185]),
array([14.191661607288673, 15.380008741034692]),
array([16.08131998959448, 15.361050594281007, 15.1620806289918, 12.142605481260127]),
array([16.53380051189933]),
array([15.285925749974037, 15.866589340311787]),
array([15.850768294546251, 15.151762982750638]),
array([15.086979817613425, 14.819950481338548, 15.870383617371791, 15.256007734192037]),
array([15.630953406051248]),
array([15.990120652030352, 15.695577393311742]),
array([15.464837975363455, 15.414209210230654, 15.032658469752583, 15.735463888635088, 14.496979519232095, 12.704316850676719, 13.816787994481734, 12.260565779793932]),
array([15.423357271944152, 15.532444491463512, 15.409841027044232]),
array([16.08537479061233, 14.186670935855675, 14.877231271032898, 14.61188701789669, 12.242103757612737, 11.800948553624618, 12.072339335327161, 13.195768966743433, 12.752108066655747]),
array([15.174255888182994, 15.750242514107368, 14.318952462335469, 14.239053030737857, 13.930273686060094, 15.51934694811125, 13.697298449832033, 12.883806942929528]),
array([13.925716375935652, 15.679980147822642, 13.840348341060038, 14.251205109567657, 13.98062234242089, 13.676572718576697]),
array([14.488833886784235, 14.072275981102717, 14.659985157954946, 11.982263244933684, 11.028729118446197]),
array([14.56777405223256, 14.211324490492856, 13.981049166120544]),
array([14.50702797186518, 13.617514654212377, 13.02716780043332]),
array([15.44264544832977, 12.884675798266525, 11.796973178502098, 12.283461760803384, 13.59422042623773, 13.417304162472998, 10.638623315777584]),
array([14.552745901034001, 14.535904834723407, 14.991310570979502, 14.033728806254802, 11.758339940846772, 13.006499752431559, 13.085898132712936, 13.073749723554458, 12.351727490066441, 9.186304006295945]),
array([13.859360501683465, 15.01414450196225, 6.648516463355727, 5.30794728147368]),
array([14.770103670926186]),
array([14.062592582859226, 14.64316660084839, 14.155020105989683, 14.231989050932354, 13.571673810829623, 12.673129179786935, 12.95470925466851, 13.105550764487244, 13.23905241648124, 10.954114113781694, 10.503794948955676, 10.426291381976418]),
array([14.459238134656573, 14.789612918897511]),
array([12.266865555074414, 12.540065413236503, 12.466524513714184, 12.970196014503296, 12.916621920252984, 10.17412045029385]),
array([14.551882620465236, 14.328345877274815]),
array([14.166864864105149, 13.819772535068383]),
array([13.646752583214786, 13.016634166721419, 13.288678505371635, 12.223593373390772, 12.72412617845832, 13.661374656292208]),
array([12.456098567906256]),
array([14.024089700066959, 13.135532506559816, 12.139893496920928, 12.593393450859775, 9.73896262235839, 8.487971147017086, 10.896278321764916]),
array([14.05366422900163, 13.999619309074676, 11.695963912355815, 13.351665138703778, 13.22078436270033, 12.886619044334601, 12.547289045552064]),
array([13.368300433605196]),
array([13.513338002049354]),
array([12.364542770666624, 8.168227033064824, 7.422538223337616, 9.551229544988974, 9.255070023226002, 6.984291655832013, 6.667985781133817]),
array([13.059486180455139, 12.218608105123979, 13.150479819887217, 13.00673954274206, 12.22247318260497, 12.164740654698804, 12.810405868594593, 13.104767008551233]),
array([13.085822687732279]),
array([13.128960245747328, 11.16278845053497, 9.363489129557925]),
array([12.038379590825075, 12.281086592632365, 11.846401702890368, 12.723444432189519, 11.636151281852]),
array([12.287099447026977, 10.484292675343038]),
array([9.571119451561813, 11.573035035360867, 5.599663104842158, 5.152947396813226]),
array([12.3766812651295]),
array([8.65211970271025]),
array([11.98567537447158]),
array([11.936650355595, 11.877362382449835]),
array([11.296473876800988]),
array([10.725338560485383]),
array([11.175553464731877]),
array([10.588201294718036, 9.620705531743166, 10.15485519570461, 9.723930661853347]),
array([8.431113226880171, 4.917342564481488, 2.0984259977341315, 0.0]),
array([7.576780826670248, 6.253986521256582]),
array([8.049005421506266, 9.55221507862095, 10.247935371512119, 10.310644415803631, 8.599765314254354, 6.291997392703543, 6.9232170154500325]),
array([9.524602463886552, 7.794543217766591, 10.057025757597224, 8.93993684820241, 8.248472593022578, 7.81181976102989, 8.71513004560761]),
array([9.243578832328474, 8.135760071342652]),
array([9.30933455949389, 7.986320858137082, 5.516751516909733])
]
d = [data_1]
names = ["45"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T43', 'T44', 'T45', 'T46', 'T47', 'T48', 'T50', 'T51', 'T55', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T85', 'T88', 'T90', 'T91', 'T92', 'T93', 'T94', 'T95', 'T96', 'T97', 'T99', 'T100', 'T101', 'T103', 'T105', 'T106', 'T107', 'T108', 'T109', 'T111', 'T114', 'T116', 'T117', 'T118', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T130', 'T131', 'T132', 'T135', 'T137', 'T138', 'T139', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T156', 'T157', 'T159', 'T160', 'T161', 'T163', 'T164', 'T165', 'T166', 'T167', 'T168', 'T169', 'T170', 'T171', 'T172', 'T173', 'T174', 'T175', 'T177', 'T178', 'T179', 'T181', 'T182', 'T186', 'T188', 'T189', 'T192', 'T194', 'T196', 'T198', 'T199', 'T204', 'T206', 'T207', 'T208', 'T209', 'T211', 'T212', 'T213', 'T217', 'T218']
def get_taxa_names(): return taxa_names