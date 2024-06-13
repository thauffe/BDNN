#!/usr/bin/env python
from numpy import *
data_1 = [
array([34.919341630506096, 30.057813679778697, 28.774548866869537, 33.628372539612414, 23.335316662697807, 27.431795298532485, 23.288745282368435, 25.545948574972183, 20.857809142088485, 21.72701135733529, 20.594151728772214, 21.726067270534696, 17.765235324672094, 20.41855101406551, 15.803785413917955, 15.669627143228256]),
array([28.354864912831804, 26.230851593391883, 25.669196355655096, 25.440705789198432]),
array([26.536246839988983]),
array([28.960374991957988, 23.144209554187647, 20.969257883845284, 21.02813059047846, 13.768830839715163, 10.122577279033072]),
array([25.932617131278867, 26.667135433380558, 23.033403393766356, 25.219332951644866, 23.87092994705648, 21.93212659710221, 20.457744882462297, 21.076569284981158, 18.737209604376197]),
array([29.058820017097556, 23.278250017587197, 25.966624456947237]),
array([29.34231834965369, 27.626597611522545, 23.949373323105693, 26.33797537641611, 24.53566123123406, 21.908816914899685]),
array([25.963939621945798, 27.13396592543747, 24.830179826425862, 20.68072031193933, 21.2983836291369, 22.60246204401709, 18.490046349149967, 18.605875027583394, 19.98678023067278]),
array([24.863139542778203, 27.621535868033714, 22.42043788916018, 21.85613091711055, 23.013460709725468, 21.44804702497071, 20.446573841725176, 19.872871895530483, 13.165480397855264, 11.680252047317182, 12.782975553550758, 13.053178790735, 11.397584409340721, 10.975832672096102, 11.578137905191774]),
array([24.068891841608234, 25.45004474542498, 25.977866645900193, 26.184740493400934, 23.48719121769388, 25.130732533217767, 27.678402507182184, 27.313773147315608, 26.90083298005069, 25.757892615897642, 27.501419676324673, 23.193705446451386, 24.585762279622127, 28.052126409627448, 25.856129340776675, 24.126566238200848, 26.21635120493038, 21.87496090412591, 20.532573893762752, 21.36118832297328, 22.589605268692218, 16.81005092745538, 16.54766132446509, 16.91159016294595, 14.996437215912797, 15.486738574362331, 14.035777840132223, 15.163247951167396, 14.674302286684357, 14.81922404216208, 13.405853458567208, 12.669147821846657, 13.180462244604978, 13.383404358235067, 13.776693313314313, 13.166699971847244, 12.97852761034567, 12.584403234287398, 12.95854359278974, 13.029307221095372, 13.604134328450066, 12.70073882716471, 12.983047458861321, 13.396827321214907, 13.067904018578457]),
array([23.72445843218478, 26.849612907760527, 22.81579259113736, 15.409725003022704, 15.280762879557013]),
array([26.201009300470947, 23.664052213939215, 26.082815240883143, 22.829814697100563, 22.26031748301984, 22.716038804615245, 22.381232285182783, 19.105713987050486, 19.474756181925216]),
array([26.38503343991798, 26.027543995083686, 23.85589365736294, 16.591868477003995, 18.541725584465706, 13.46159520270869, 13.029451757633074, 11.791307590531828, 12.676544578918568, 11.151613239385801, 10.505960016904872, 7.502577030310431, 7.858153647682992, 9.40587522363909, 9.063115406511685, 8.152709716054817]),
array([27.014775865474668, 24.253879297088513, 25.29905709584432, 23.776060785562493, 23.392976923768753, 24.022546689862434, 21.11009929854049, 21.029910324452338, 22.580463386445146, 20.955253242356324, 22.772154951014475, 22.891978487287382, 21.170663468330325, 21.30464006899784]),
array([25.299836992637744, 26.523297342786183, 25.673617376054782, 25.14507978320441, 26.879404043130982, 24.62038278940306, 26.282611686173464, 21.170698820607456, 22.902821604120398, 22.519941467669568, 20.83089047214993, 22.550920773202883, 22.74371290314568, 20.505462999429405, 17.43764997761837, 19.71011347620187, 18.068028717721482, 16.47740043761668, 16.088230581100852, 17.93427964355979, 17.886550766284355, 16.19965322092679, 18.80754114491368, 14.62441477654999, 14.735175252534175, 15.235300550582158, 15.602190087638022]),
array([25.310519072311898, 23.4621114119889, 25.974576945237633, 22.06520228403937, 22.71269904416144, 20.74385752974624, 19.546044161292773]),
array([23.352590049531777]),
array([14.18461831704844, 11.084114675511332, 7.795877144801519]),
array([24.33940038636894, 24.950518669110348]),
array([24.08869702524026, 23.695958632863093, 24.57445422996314, 25.466690706683327, 23.22400562163987, 21.296629372742096, 20.867347716709038, 21.18166017829866, 17.428515533786303, 20.338484167645525]),
array([17.62068352749559]),
array([24.605992887197456, 22.741213051124845]),
array([24.11244895882406, 24.207402718015363, 24.408734305399545, 23.405981148695517, 24.220233547711434, 23.800832736707765, 23.40340786977742, 22.611260800333362, 22.35317938735713]),
array([22.60752480375523, 21.117652137492883]),
array([22.76964022029298, 22.817686294166595, 22.89867334559086, 19.323679275106723, 20.263823693783355, 20.26645315597935]),
array([23.900283026583832, 22.85759182637397, 20.80671386315092, 21.579286479441528, 21.615870511166694, 22.785052906454098, 20.65966652992462, 20.731627594046845, 18.120818348342773]),
array([24.344122429697954, 23.772976126038596, 23.946211387211175, 23.50738615761298, 23.06787974295039, 24.355182393430333, 24.603042024718764, 24.05366373730983, 21.717579859781775, 20.70291871377615, 22.34235142325546, 22.969777300306227, 19.349036813280204]),
array([24.953193689106012, 24.893730964857497, 24.719786566313037]),
array([23.39044639710464, 24.58434813407094, 23.75836667462359, 23.992587715299514, 23.69550303620942, 23.269139892378718, 21.947967538394273, 22.50957861338028, 16.514155200515475, 16.781241297646496, 17.241430072722924, 20.274306220943007, 20.001353204173263, 13.834245412187435, 13.721550219236326, 13.272909996073139, 12.999767853258703, 13.52880458479255, 13.524533164559223]),
array([23.1881706972086, 24.51648710374885, 24.479457970701713, 24.382863855379465, 21.066397520286486, 22.69307486102347, 22.0231504590589, 20.269871239839603]),
array([24.169141139338482, 21.61792207546575, 16.278270692186535, 19.505972060229123, 15.517230193782353]),
array([24.544335509581757]),
array([21.572974934308437, 22.99985604157592, 22.631706539818623, 16.804738046812467, 16.86699135422531, 18.334739921999137, 20.20395604876623, 15.703682588717038, 13.679956892973838, 13.670355816366389, 13.418548926834422]),
array([23.72546570322178]),
array([24.388541858704077, 24.192932985319185]),
array([23.12883536961103, 23.51764907551422, 23.294945356574534, 24.070393844118154, 22.160626047418866, 22.90501459366511, 22.17629349343069, 20.937404434908817, 20.68160377911551, 22.296947523913087, 22.444239896662268, 16.222906054986098, 16.885556581014885, 20.291283723274915, 13.856211197439595, 14.44351017495524, 14.5939406930545, 13.22345189973744, 12.333849837675336, 12.964313247898234, 12.971367300314128]),
array([19.62720876237051, 16.133627349086225]),
array([22.37664860662359, 22.988212457393654, 22.809833938254005, 20.796246999703772, 18.713378348085673, 19.628111605197002]),
array([23.741125945818123, 22.879237556059834, 20.4632665018098, 18.644975687301017, 18.923442396452245]),
array([22.62280826910062]),
array([20.538833843583394, 21.621870617171002]),
array([23.64167586728126, 23.50862844513921, 20.45775657754805, 20.7741992540742, 22.632906176017386, 21.45918353195744, 15.98613990730314, 18.86006622993377, 18.153139516357314, 16.16899957139301, 16.88827008843902, 19.746235853808876, 17.684526066663395, 14.355328359808823, 14.583921844607655, 15.425908596947263, 15.278735951126913, 14.934748548302503, 13.778696825922635, 13.660811031417072, 13.720798474462228, 12.926587361461653, 13.400685529225395, 12.942175076956602]),
array([20.55584003832334, 22.202849741703158, 21.639774635670403, 19.917347549933094, 19.536542518160363, 16.9490653110326, 19.18322867592992, 15.29409262568607, 15.846516132449251, 14.81558070247452, 12.875704753508744, 12.628799645592007, 12.498482381025926, 13.708197415336052]),
array([23.69162354792175, 22.55395680215365, 20.668470401012904, 16.59207227571889, 19.8860050428328, 19.612316618706004, 18.114677163133265]),
array([21.262154727734472]),
array([21.69681517769129, 22.481173775367253, 18.77198346627789, 16.039120355209, 17.800774340025846, 19.96469557337221, 20.177508122199736, 19.008401014516572]),
array([21.9271453432257, 21.411155894112017, 22.134058422672318, 22.14030601804946, 16.707025378311002, 20.23479573487045, 19.361971061585677, 17.026121734267697, 19.960994949123943]),
array([23.07138127695683, 21.53432159561781, 19.680771364427496, 19.140722201414597]),
array([20.670812949177375, 18.239374650588886]),
array([23.230861939500127, 22.493937226417582, 21.292048921411226, 17.261242864497994, 20.13888271284293, 16.502508536521503, 16.00232542278548, 18.921553795714377, 12.535137124919, 12.960982394231578, 12.770850912470094, 13.646295153781816, 12.617992175651253, 13.645816016797742, 11.111569383575086, 10.881535287732776, 11.476299732592034, 10.779301326289202, 10.78395749056804]),
array([21.89711417498481, 15.479136136796287, 15.045706256232524, 15.284531133247965]),
array([20.723876245171834, 22.17203067847362, 20.73002804094931, 21.69278769168321, 21.99866390194575, 17.97173004170733, 19.809384899865012, 20.109878729782185, 18.97645536365827, 15.066155387585722, 14.405766677654494, 14.361748212148978, 14.738825486485608, 14.410054542249886, 15.947019601740768, 15.503543444663544, 14.174969718315417]),
array([22.41075320221601, 22.907085549048308, 17.497321593858032, 15.292459285567679, 12.886695413164281, 12.10390822646032, 13.038822274843069, 12.733397812638199, 13.195160908648223, 11.975207132557943, 12.706672944069322]),
array([20.68718640100443, 20.970996895612586, 21.965500684519725, 21.126503596861898, 20.203525964430803, 18.4050909727756, 16.770135583862167, 15.926010418620137, 14.024393828622193, 14.819015621650326, 13.713975322850352, 11.683837486742508, 12.165569938782852, 10.138651868900338, 8.842461411379368, 9.229339821419877, 11.319449131730227, 9.738427720167515, 8.797018841300634, 9.903740415817248, 9.69109830519243, 11.035371044208748]),
array([20.770651559951855, 20.611775465545332, 20.675590732246068, 21.8164939023605, 22.057460645365893, 22.704053272123303, 20.932015956188206, 21.1975267999614, 22.050387226014315, 21.88340572837074, 17.48936236564094, 19.181378956477378, 17.405008219090142, 19.317972214858575, 18.98025183502307, 14.59138618036732, 14.25116964196755, 15.642150512363253, 13.62269621466523, 13.552094178860008, 13.54295732570683, 13.120336655432396, 12.84274507118157, 13.626962845163886, 13.383445548923989]),
array([20.461976715897965, 20.627319319613044, 22.73424087936558, 17.184626284646388, 16.04230550881215, 15.285892882413606, 15.676363567571322, 15.809984318381739, 15.285798867391227]),
array([22.552549717863435, 21.61445302572492]),
array([19.910906159553175]),
array([20.831324229279257, 22.490063180053923, 16.356822249637585, 14.669467338152375, 15.580423983148226]),
array([22.459172486467427, 16.18770116824341, 17.78819108218304, 16.82344946838589, 15.729263910668966, 15.155182518690628, 15.120111741089358, 12.910906918230477, 13.051230529293415, 13.107014338581239, 12.948201464775705, 12.043826187214265, 12.510453740758079, 7.653233585714794, 9.252426830342388, 10.372234283069995, 7.9419070882434575, 11.160358121322759, 7.785012476043614, 7.276362869868016, 8.077847989642029, 9.732967527247744, 7.277337284012732, 6.670130704050298, 7.049856321066091, 6.043337495062754, 6.630545766707767, 6.366996712608854, 4.126662414778282, 4.42786859956936, 4.474967450737366, 4.087465236396426]),
array([22.158713420688958, 18.462963671499118, 18.714167400015445, 20.341393260864145, 19.100615571731428, 19.92729945196102, 15.870027466583867]),
array([21.498288467366155, 20.52748069250765, 18.717581473965843, 17.907354535728654, 15.812725943831044, 14.66879907197663, 15.720722258393552, 13.651843525819775, 13.283474759288556, 13.249796418893663]),
array([22.125255135048928, 21.2301447595439, 20.95965931259383, 15.657373649742507]),
array([21.887689980116942, 20.132530209404408, 19.514370832893736, 19.48242028884963]),
array([20.556823058481573, 16.935715677776223, 16.304149563401634, 15.925845226626631, 15.250606699141102, 11.70828980915849, 12.840656643255588, 11.4678818407296]),
array([20.818168379468524, 21.139793783783247, 21.65224974507447, 18.897867646189002, 16.06468718243415, 19.303835211052043, 17.712457289593193, 16.69437825349551, 15.61751517422346, 14.599839499205038, 14.652522761590108, 13.856735929966742, 13.849169305046992, 13.768844640833251]),
array([21.02798971366825, 21.126716011779312, 17.077042551831745, 18.37697113990319, 18.936042779561056, 20.169908623826455, 16.583595207996748, 15.097477150401257, 15.369163139857093, 15.134692532273958, 14.207171311315326, 13.617839563530707, 13.524781434175377, 13.379505493071957]),
array([20.540084376402206, 20.4730425764373, 16.847001201854454, 16.27590666922802, 19.94830794779172, 16.939220676558346, 19.96850160640558, 20.326140752488186, 20.112735615972213, 20.274634249866693, 13.920148523365024, 15.098998877108235, 14.816822952949014, 13.696605103921936, 13.631871999985114]),
array([21.34706657194033]),
array([21.03656967473997, 21.46972585683253, 19.14049923149496, 17.452945713266736, 19.08334941198879, 17.200584857763715, 16.532683957743316, 14.657683802266183, 14.881817170457985, 13.539459576074629, 12.485315346918462, 12.294726171294645, 12.85589494043921, 12.142749318454683, 12.46921022649766]),
array([20.513657992228307, 20.755469653403715, 21.308731198418236, 17.847558726467845, 17.617497573186, 16.5284476602365, 16.125400995166096]),
array([20.46554665627673, 18.821556940362196, 18.290734690376677, 20.15370729708844, 20.340698736323223, 16.785160511470846, 16.681294480670076, 20.15594063760314, 18.011250947874732, 19.456684129308844, 20.03415156057501, 19.286419810868125, 14.120832699494292, 13.942556802475957, 14.194156340040017, 15.024883957945416, 15.576240243989572, 15.584568723772973, 14.538026772037039, 15.265564672089713, 15.290826594977723, 14.740518424080452, 12.092246064287286, 13.005565632738083, 12.916195168138094, 13.436788272995587, 12.12077969513202, 13.013880500431133, 12.833561506990776, 13.298713909494918, 13.285086905928809, 12.950417677583843, 13.383624922836438, 13.588192799538916, 12.178637214891154, 13.777195235140791, 12.006499163145763]),
array([11.841296446301925, 13.532340056255862, 12.73278455947939]),
array([21.13443766247377]),
array([17.69129908744796, 20.011138691213652, 18.347265371259425, 15.968769338171914, 12.513288976942835, 13.084634665208508, 12.687157548648344, 12.30690520455023, 13.423571269214571, 11.408864225125289, 11.369514270279243, 11.422589142369915]),
array([18.034809859775244, 17.748000230699944, 17.616518782870962, 17.936538993107295]),
array([18.364732480154686, 18.250000368901876, 16.200048035050948, 14.631072795320707]),
array([20.721556314231165, 20.86821492553613, 17.541950591887492, 16.925046568942452, 18.249126567538664, 18.63058513528093]),
array([18.218036056737542]),
array([15.694063365883366]),
array([20.54943679064794, 16.06100892939119, 16.26851322825257]),
array([15.335564949680734]),
array([20.444337039718377]),
array([18.764150982264688, 16.211317634634646, 18.441399027418196, 16.850594450213226, 17.53856484723276, 13.32067032825382, 13.15224061857788, 11.911632821889446, 13.28016029749473, 12.788919651600304, 13.083952805687094, 8.66813444687266, 8.098471108312507, 10.98822391680688, 10.45218051326266, 9.291494101389006, 7.879596538398669]),
array([20.51511818703468, 19.51327717281092, 20.268516002465642, 20.30835072794256, 18.709591143155592]),
array([18.542536016162003]),
array([17.844514856724835, 18.58291747487888, 19.314678083595442]),
array([13.870620981971141, 13.814353348704634, 12.644678278265287]),
array([19.117951513031816, 18.026385262949358, 17.973734981152518, 17.002786130354576, 16.653533796355262, 17.059521122050516, 14.73927124006469, 14.429411124422495, 14.27842581111906, 15.023882039417815, 14.857271245656635, 15.952110010869164, 13.812376040193094]),
array([19.59759405081054, 19.03749108341176, 14.035514441289171, 15.58470617651886, 15.058313559865, 14.698579731744521, 13.60669823095566]),
array([15.674863141913193]),
array([16.97899008840368, 20.018721303210896]),
array([19.430404777645606, 19.268903040602037]),
array([16.38166303584891, 16.851315086554024, 17.50497653294219, 17.905621624452916, 14.889069954822443, 15.070242558669442, 14.989121202123346, 14.121898001975813, 14.70821778810914, 14.395486885888028, 15.49878645075742, 13.369194189698307, 13.459902185417018, 13.618386379277723, 13.66296020630082]),
array([17.875149384696186, 18.56265451069074, 17.441241456133902, 12.564767991746004, 13.664510620855165, 13.719865264493322, 13.377349150045339, 12.851282678602512, 11.601558404842953, 10.59800977240215]),
array([19.083877914572916, 16.45901258239159, 19.019026615475603, 14.63105787077125, 12.551932703262704, 13.637408479186636, 12.471384545200788, 12.355217338324321]),
array([17.559348068083914, 19.194162956581874, 12.346147156025717, 11.838734658496172, 12.540998390340702, 12.89878958448249, 13.354624915726072, 13.563436728983579, 13.663498056665619, 13.352759084455185, 12.632256752587132, 12.307648462139882, 12.065639471373505, 11.61496162996077]),
array([17.63142022235472, 16.516759853894683, 14.455546987273648, 13.482834359381325, 12.889400580721967, 10.064481894463952, 11.584369098749352, 9.951782271984548]),
array([18.19425651274263]),
array([17.289884590012715, 18.617065100925593]),
array([16.05959081335835, 17.002511685315064, 18.65092410497848, 17.01418350265248, 16.022631129326534, 19.129002256739422, 18.13913741788402]),
array([18.565008461968862]),
array([16.13191453653461]),
array([17.931728616827638, 17.77220019841431, 18.762911595045544]),
array([18.65284778310992, 14.676784821192578, 14.880719428286218]),
array([17.613881190041344, 19.342321298486876, 17.476342733233725, 18.937915558530303]),
array([18.46404448541064, 19.152190133864956]),
array([15.76018313947839, 15.459535549338655, 12.704609450809906, 12.193062140506703, 12.390983264522314, 10.549496102953372, 7.592122989898138, 8.751103326826504, 9.477343496080358, 11.02083173133567, 9.976171502021392, 10.777105900083559, 9.066429059166907]),
array([18.553375969048826, 19.517574261378922, 17.600064423009005, 14.976122754322109, 12.762698136053034, 13.243268525972058, 11.792611112545158, 13.72675463072984, 11.398262023696354]),
array([17.56972376530379, 18.266585767265056, 17.687587535377055, 18.5341014345996]),
array([17.980694264862993, 17.92750516942564, 15.300134963944897, 14.229416935282623, 14.946206072166824, 15.906521919525668, 15.337620412455564, 13.062160515160926, 13.23760808563655, 13.029645976216903, 13.554671657808651, 13.529904622404738, 13.257105143770058]),
array([15.063380619334302, 15.822305506267652, 13.332768006181457]),
array([17.886966071837577]),
array([18.856002167051994]),
array([18.02980150168931, 18.75847200075012, 14.189467865132567, 15.651517973938061, 15.657005014933233, 14.33015432394924, 14.948987904167206, 15.504449442017751, 15.890642367389583, 12.851373416285876, 13.537467890903436, 13.409051376287378, 13.292137284578954]),
array([18.942290838491825, 17.233920427129163, 18.715584348419533, 16.68439191112718, 15.455539131508212, 14.537694830394328, 14.259294996304417, 14.583181878927196, 14.540130067732154, 13.487503516890136, 13.6752193499795, 13.318804628403324]),
array([17.97368441971503, 16.462031235520907, 16.075362946229834, 16.607744109853062, 16.18945421493962]),
array([15.454082504405537, 15.881511021977236, 15.630059072463224, 14.650526790961337, 11.702971567106196, 12.258134810674727, 13.396997036921883, 13.45354267375404, 12.683181219569867, 12.036083904048407, 11.701784713652627, 11.127075981068176]),
array([18.39840233191583, 18.732610907752573, 17.236286088096794, 14.947326788799382, 14.508777528608526, 15.232407137279537]),
array([17.57511822957793, 16.680366328939538, 16.819411329626686, 14.576678707129089, 14.25890302512627, 15.341147583295536, 15.836639212657065, 15.937058326743722, 14.87093715420575, 12.293645193769814, 13.035446588416775, 12.24511641815289, 12.684911282274912, 12.575454423516137, 11.757415673600146, 11.738133042969604, 10.73080333164609, 7.474987519300041, 7.471134093071787, 8.364519469161934, 8.684900897285436, 7.711710789768546, 10.925298623334099, 8.817491041902166, 11.480815202364568, 11.451970786220246, 9.913520817231984, 6.96090439529645]),
array([18.46144708759342, 18.685531047126666, 15.577791499207848, 14.648472635121342, 14.511461284629126, 12.582750596212955, 11.048964953338034, 11.016031516383867, 8.672489850704595, 9.979697378380582, 9.956303879508104, 8.956671793056998, 9.594376250273458, 11.110841248891163, 10.216880455965551, 10.429429388555713]),
array([17.692450539150038, 17.643669989676894, 17.66446378204597, 17.798415629750952]),
array([17.982473866125428, 17.467590324464034, 18.251026671076414, 13.790677621074984, 13.546227991240926, 13.190042831197129]),
array([11.029363473531996, 11.30727888490131]),
array([17.856614432618443, 17.87277698790829, 11.642417185217663, 12.827423536177537, 13.092286329580348, 12.743118709299166, 8.852279939312005, 8.602333829901477, 8.9482923312895, 10.998274399871702, 11.255447909234482, 10.319089919902348, 9.92323782962112, 9.327499710545196]),
array([14.58848659808134, 13.548605825258262]),
array([18.073060206008417, 17.54273257171123, 13.966096123639192, 15.907510998699895, 12.311898237907128, 13.313878458623932, 13.689689735677257, 10.167843763634796, 10.977744739939173, 11.04725573375539]),
array([17.12668379358621, 17.15913520694867, 16.734321399178384, 15.281431958960578]),
array([16.049103785536424, 16.189612032507533, 15.705481545648755]),
array([17.26245042355613, 17.93773948003542]),
array([12.47508138789284]),
array([18.188033464360537, 18.28392418992159]),
array([16.43995576779967, 17.675229042334436, 15.987786923785352, 16.425269918847707, 17.417723950086017, 15.884538977099492, 15.88073256558439, 14.788316124128084]),
array([17.57058766679085, 17.051964576198024]),
array([16.393339080786305]),
array([15.591046276489914]),
array([16.311973788566732, 16.242993935779904, 14.657594203338217, 13.84117027052134, 15.571079339468579, 14.049489017907076, 15.480567250168855, 13.5080261328996, 12.596116411186806, 13.740849439688432, 12.51983895802235]),
array([17.111658328891394, 14.651828024629014, 14.188146613414453, 12.994935315910759, 11.932737299778882, 11.873618285007364, 12.610761065028512, 13.402342782151733, 13.516442429737996, 12.83986791008991, 12.358529866542655, 12.474651252768481, 13.45734242408103, 12.037985321855755, 10.978055286483567, 11.220204566809866, 10.848225766797123, 11.606504594592453]),
array([17.45990087228898, 17.675767908802893, 15.287602040506446, 15.679373881600235, 14.934713972393649, 14.492458246526212, 14.69686216173793, 15.333677535430565, 15.39194013475635, 13.52242969341196, 13.327063852780046, 13.469772497991642, 13.477735535539942]),
array([17.012265590283302, 16.53594444550811, 15.350857713405835, 14.250282128429596, 12.734430394986083, 11.797658116300514, 11.926918471367944, 12.580220818101058, 12.993057308894215, 12.177296491204771, 12.89180726759562, 12.886834897587713, 11.413185661515932, 11.306452817669218, 11.210248765116102, 11.293780216931387, 11.479997887685936]),
array([15.210851364644688, 15.047547859488422, 15.758084933685021, 14.482087785748966, 15.23636145856615, 15.853273749065734, 15.943387862826347, 15.721067580433226, 13.638217622306382, 13.73948899792268, 13.22282779122403, 13.30659479338225, 13.712432982060115, 13.559880036814668, 10.816383357652796, 9.360735922821721, 7.524674492103903, 8.059531946973035, 10.009887285561938, 8.108476819635857, 10.502303982679525, 8.118892023747925, 7.642594591392225, 7.626044755153121, 7.285697832468922, 7.986621363784539, 8.681818833642588, 7.834001912120849, 10.003393329509905, 9.929545294263987]),
array([16.414374104818016, 12.501659569155185, 12.322880918940216, 11.644086019911201, 11.550739804704715]),
array([16.31403652963034, 15.269062134467125, 13.837553671981983, 13.143365719013103]),
array([13.315577689064792, 12.014056485127771]),
array([16.915808626624084, 14.39667500399764, 14.38617154820744, 13.896394555464417, 13.420169405067359, 12.592294411318674, 11.778280041152748, 13.44508795991995, 12.686755380236182, 12.05985137822336, 11.246955767797894, 11.285338343287586, 11.456037211535524]),
array([17.141711565635287]),
array([17.224275733413513, 14.86877680899799, 14.614641415023081, 12.522097272307127, 13.240715534269283, 13.043102178870349, 12.28273862765564, 11.288512422243436, 9.786957831424761, 9.994245934344223, 10.055032567325881, 10.19945237556265, 11.621842258841914, 9.774548214873084, 10.1240489656459, 9.521570537093973, 10.445065308727774]),
array([16.240866503976243, 13.431528307234705, 12.653002922415913, 12.908335209322543, 13.057183419016877, 9.640787586226278, 10.644440667452413, 11.102154748359444, 11.141405571511443, 9.352364903628423, 11.408691275540459, 10.16329358412152, 10.941936231147087, 10.894061731128142]),
array([14.485468313357634, 15.249967538031985, 15.706449782067002, 14.576312043271205, 15.378767667293173]),
array([15.106268935235557]),
array([16.997062222904827, 14.652844317737786, 14.932442673183862]),
array([11.76492868377955, 13.081231585157985, 9.996357435296165, 9.9863698571677, 10.78216488558414]),
array([14.848471730176716, 15.21863135783971, 15.655687455147332, 14.787070294851464, 15.638409535976216, 14.659818792957624]),
array([16.526602867174404, 15.915668907545923, 15.644679793474046, 14.873800915021587, 13.373226413806849, 11.821565025458236, 12.277639966862221, 11.856177291859797, 12.237324935672428, 12.466224594175468, 8.550368553903432, 10.01276810791896, 10.944681608238412, 10.324212079718322, 8.828102729970013, 9.567285890132034, 10.275379959487212]),
array([14.013554417527285, 14.30074087414048, 13.876304468507547, 13.462835816172788, 12.405770279738327, 13.455438256803554, 12.23529802674639, 12.70088481014799, 11.288708563268601, 10.210277785942596]),
array([13.255214789639057, 11.594576800512717, 1.6005878857363347, 0.0]),
array([16.29787193773266, 14.201370348990327, 14.883464376243994, 14.680649377499345, 12.524591221845132, 13.348749597963801, 13.371658701724652, 13.151129227753321, 12.73310572282917, 12.437180752815008, 10.872034648539307]),
array([15.172361087975833]),
array([15.452471336009557, 13.813574708774118, 13.466675020275158, 13.767603629134424, 12.336319929731197]),
array([15.235012423322253, 12.962126228085845, 13.094672482009612, 13.603674972314098]),
array([15.113941224399824, 13.813836347138801, 13.678313988568602]),
array([15.698149300223355, 14.74886309820967, 13.188299348281168, 13.02275739010107, 13.123766048003764, 13.34965311371832, 13.288450193633581, 12.837693587765013, 10.077482953891334, 10.985433657014422, 10.892032300831799, 9.674344136600647, 11.21169726084345, 11.25975082653574, 9.704088220147778, 10.12171981247557, 10.202614294342288]),
array([14.119357386701362, 15.638366941879069, 15.318091498146153, 13.459434263092346, 13.666871188431712, 13.408083210450231]),
array([14.07867518977092, 14.926790496435911, 13.605678969089018, 13.489449370049993, 13.633839212256309, 13.578696429674382, 13.662653098248501, 13.532425296180094]),
array([15.015946777115163, 15.260439152013506]),
array([14.007565916280027, 14.864350177252943, 15.36849232294249, 13.176396244487233, 12.809767552267054, 11.730176880629568]),
array([15.520598021762083, 15.393100653701483]),
array([14.534896011633014, 13.883049914231789, 15.158399171734706, 13.76798391304647, 12.855857047201152, 13.316190317148221, 12.76530723860974, 12.863816178191865, 11.969657661259774]),
array([14.821407484926478, 13.410033174676968]),
array([15.214000643572342]),
array([14.37369048719202, 14.484528469375995, 14.760735851152518]),
array([14.346874035714238, 12.820336072740076, 12.759521513004666]),
array([14.466178953623352, 14.803534078776849, 14.676801054946154, 14.126097546107143, 13.773548838530518, 13.158377932074888]),
array([14.32400838788459, 12.317137714131107, 12.139360291293398, 11.072780252907593]),
array([13.373827474583365, 13.78896928391131, 12.596142459799399, 11.966371454712974, 13.066865879011063, 12.439415212308933, 8.949882102709608, 10.703328568132216, 9.637040662555254, 11.421993639176762, 10.199853881322696, 9.154946176028744, 9.621121376942536, 8.359596861900467]),
array([14.044192675667215, 12.91085281269405, 10.26631867354472, 11.620604413502823, 8.024879476668765, 9.139868415210398, 9.60669210486306, 6.030278206785912, 4.367654303122556, 5.2176821306545635, 4.903115051499277, 5.263722688778663]),
array([13.576728529955128]),
array([13.254083134650799]),
array([13.045144940430905, 12.664262943299892, 13.279550855024441, 12.683496728279565, 10.123745681827259, 9.853768597753966, 9.974280680581554, 10.012204109404891]),
array([11.711678333286828, 13.468694480714891, 12.6312793878968, 11.933743430655275, 13.465129975355115, 11.947737995044472, 13.009639029812265, 12.073333544302276, 12.841514954327447, 13.065227968847825, 11.749674458858575, 11.872134562189371, 13.332912333720463, 13.49142605455228, 13.160769085555165, 12.33497162582822, 11.202887475238054, 11.191542247885211, 10.380634820848028, 11.471392316799376, 11.041926066264173, 11.40124600062977]),
array([12.355119021079684]),
array([13.466908064505105]),
array([13.441427858544621, 12.441321238635343, 10.17450358229526, 10.016242214660535]),
array([12.509446771896418, 11.956940145263992]),
array([11.808310372835013, 12.933509866919433, 11.630306988513052, 10.473653493100215, 8.756580022414962, 11.491923158581045, 9.198381256890604, 10.364469829240255, 10.41644256212004, 9.089788976213597, 10.088072167166859, 9.20301367901282, 10.03273407512169, 11.3207446532334]),
array([12.33222631951444, 11.972148074974633, 12.343530480726969, 11.756468920299465, 11.756811289379824, 11.990416166229316, 12.859824163333657, 12.57868047939383, 12.360742280756284, 11.850234381218849, 12.209105886376086, 11.436548576795238, 11.605388139051097]),
array([12.357381169642371, 12.363330472477786, 12.335770416830462, 12.445471979104594]),
array([11.958710722184358, 12.182379691737914, 12.097223311930946, 11.865102039021332, 8.997378432227336, 11.380853613528997, 11.600439926057094, 9.825841617532138, 8.649760955948746, 8.920088350030609, 10.539071792861801, 9.419949852847125, 9.940057807551419, 7.801610960532731, 9.065011380869848, 10.323725263404022]),
array([11.473092575798278, 10.36269882353402]),
array([11.71884531364573, 11.874791297133985, 8.259607148742209, 11.092379244068292, 11.230609851923182, 9.475808218412372, 8.647625396417585, 7.966681015307363, 8.720530951496759, 7.601183602532331, 9.75780129291991, 9.228268476468374, 11.527817008952244, 8.588366762732726, 10.366326284883078, 9.14672692647827, 9.869327322531518, 8.113319774836546, 10.652712668614814, 11.139780896047933, 10.525340112622532, 10.972026020221115, 11.241107575668627, 8.473779588300893, 6.480564914112321, 6.503477347310275, 7.100521420852458, 7.166785525172807, 6.511550723524412, 5.743807589159928, 7.184346249517288, 5.902928158406887, 5.559748464856725, 5.305870203908796, 5.268982970175061, 5.2098898813602705, 4.764901541945234, 4.895108682285169, 4.41290744653584, 5.211960307964345, 4.7127104009076435, 5.300560079915862]),
array([11.967932601484344, 12.118087028009507, 10.134428743894825, 9.938168444901425, 9.608282888576879, 8.927556447430764, 10.83226403355713, 10.82883820075241]),
array([11.067255683945826, 11.55287828268922, 10.64610155275761, 11.292571135603364]),
array([11.712939878211847, 9.113663648480708, 9.703761567492409, 11.1502460588603, 10.630057327015283, 10.662936778132437, 9.046691903662095, 10.40315362548209, 10.847192658851984]),
array([7.77614137186341, 8.657047412768016, 7.844640985595725, 10.91316324085849, 8.124302200173446, 9.214895313546387, 8.722257398578027, 8.908173651177126, 9.723986779859311]),
array([9.416159896432504, 9.184236898176028, 8.28600808826696, 9.548954133361008, 7.982622505183173, 8.415943013129498, 10.656431967633395, 7.75588006610083, 10.854396145698953, 9.04369198368757, 8.994557146908306, 8.908416645390782, 8.055532225823418, 9.217863474947023, 11.033812292622686]),
array([7.63729644666921, 7.387848290878983, 8.205229383309906, 9.216234483743317, 9.244546868618146, 8.736047108652642, 8.299379590842364, 8.941412466265296, 8.296259182384919]),
array([9.074676639369105, 8.745846782626849, 9.357448476752786, 8.293015046547094, 7.9564061926479, 9.311948526813113]),
array([6.263748456151823]),
array([6.4819730674880445]),
array([6.454044053377558]),
array([5.582999731915475, 5.446190466967958, 5.292802169274382, 5.158570933427145, 4.972236951671735, 5.2590755778343805, 5.162523416650888]),
array([2.0765594637846734, 2.4016297474537835, 2.1690899218453246])
]
d = [data_1]
names = ["91"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T16', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'T36', 'T37', 'T38', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T46', 'T47', 'T48', 'T50', 'T51', 'T52', 'T53', 'T54', 'T57', 'T58', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T83', 'T84', 'T86', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T93', 'T94', 'T95', 'T96', 'T98', 'T99', 'T100', 'T102', 'T103', 'T104', 'T106', 'T107', 'T108', 'T110', 'T111', 'T112', 'T113', 'T114', 'T115', 'T116', 'T117', 'T118', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T128', 'T129', 'T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T137', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'T144', 'T146', 'T148', 'T149', 'T151', 'T152', 'T153', 'T154', 'T156', 'T157', 'T159', 'T160', 'T162', 'T165', 'T166', 'T167', 'T170', 'T171', 'T172', 'T173', 'T174', 'T176', 'T177', 'T179', 'T180', 'T182', 'T183', 'T185', 'T186', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T195', 'T197', 'T198', 'T199', 'T202', 'T203', 'T204', 'T205', 'T207', 'T208', 'T209', 'T210', 'T211', 'T213', 'T214', 'T216', 'T219', 'T220', 'T221', 'T222', 'T223', 'T224', 'T225', 'T226', 'T227', 'T230', 'T231', 'T232', 'T233', 'T235', 'T236', 'T237', 'T238', 'T239', 'T240', 'T241', 'T242', 'T243', 'T244', 'T246']
def get_taxa_names(): return taxa_names