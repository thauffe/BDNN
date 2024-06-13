#!/usr/bin/env python
from numpy import *
data_1 = [
array([32.72064040991912, 31.58310015839438, 33.57397184447023, 33.42708092733712, 34.562083991601064, 34.962402451611304, 34.262211353743375, 34.37204118139456]),
array([31.917201344858817, 32.16716637468117, 30.135996695289535, 32.001639239115725, 30.251200615976067, 31.937092321282307, 28.309168254517267, 28.585748364071133, 31.846245850393473, 28.64930694988013, 29.221322455812473, 28.245584385678647, 30.623278320747172, 31.979005124589115, 32.9490176312657, 32.59332178552253, 29.74581480078747, 32.45508484361912, 28.545459110939838, 32.81711652982982, 29.78538660041632, 30.953295581925136, 32.84598464186712, 32.4385061413936, 32.50376110634747, 30.295341906593986, 28.549372753868827, 32.303001919463, 29.061559145781693, 28.395364758660833, 30.22336021907479, 32.73871535939133, 30.745555792516527, 32.93330813764352, 29.558940630737055, 32.920418685165004, 32.2550265437097, 32.3701167460604, 28.75177021909597, 31.797452302822077, 29.34647231158027, 32.368246823001044, 31.93511927434826, 28.9061868082119, 28.990523424746446, 30.267198594933035, 31.837666968269566, 29.43335912298591, 29.529098774995305, 30.551424858610677, 32.497336352884, 32.32725908121281, 29.37241152518493, 28.796921767144127, 32.37705851749647, 28.82899005642064, 29.284819423099147, 31.673224736411107, 30.940622954903077, 31.881261369087568, 31.463301619223984, 31.03569274275973, 32.528888738797406, 32.9289113360739, 29.879261504890685, 32.73418483711356, 31.29014999694312, 28.6134341369913, 28.47725624413169, 32.264578567948675, 31.688074636095557, 31.223382530365413, 28.33574989217887, 29.294175183070372, 30.82694108187316, 28.831987298149212, 29.028549177735666, 30.734548992802527, 26.980690704636306, 26.348117722352065, 26.697930078592577, 25.323591751693975, 27.986676975662558, 25.4781276645211, 25.190828927876456, 23.24131404448405, 27.69074075532819, 27.50457514287908, 26.0795219665203, 25.515938443300016, 27.541382684759693, 23.56024165740086, 26.533787770797385, 23.869658470316864, 25.010268489766936, 23.663417227623874, 26.05793621209831, 25.07719261114981, 24.670304568801782, 23.20352382279616, 24.33322253046252, 27.85945114111161, 24.23192643059567, 25.3388340405761, 27.169860701507286, 23.639655991940227, 21.170835749358176, 21.619003901350688, 22.620625302331128, 21.36565515364846, 21.40468289910383, 22.67922754038025, 22.075190929590022, 21.116963896123337]),
array([31.0368730577218, 30.89883577182111, 30.34828477793326, 30.086111335600926, 29.25306419697269, 28.37557364259878, 29.899170147765755, 31.7950319985609, 30.61946701909766, 28.93365214983118, 28.196524523953954, 31.352305842631065, 32.04637224952993, 28.638438573515568, 30.304424062915, 29.238734224028946, 23.113849908900306, 24.21099523884135, 24.585576544446525, 20.93228030292583, 19.45435117617685, 14.95064286119999, 15.408054723810157, 12.39393194325987, 12.130681101049074, 13.080935550334686, 12.44712161095017, 6.907341241942641, 5.713840759262925, 4.446438249163666]),
array([29.989941161213178, 28.325652605612117, 29.22509102619403, 29.40772981034312, 29.29032650617348, 27.26434963911064, 26.258828204037755, 23.354138007536193, 20.466406162137503, 20.496128958070763, 19.79928554480069, 14.962119658175117, 11.621213086962806, 10.99106439402952, 5.93449076658084]),
array([26.753099214585326, 25.655738487535665, 25.35603774624568, 23.311911989602145, 24.932124096889456, 25.006904707113552, 25.78815535629555, 23.409345458200445, 26.752748463504314, 27.205472510551488, 24.95264892420338, 26.11083398532862, 25.286334417138484, 26.539077352578094, 27.068111196605287, 24.62087266133101, 22.649191283002875, 21.767774100523386, 21.74705962644453, 19.79503218193457, 18.970306847294164]),
array([24.908227542967825, 23.778555696182465, 26.3165784750565, 25.349867933305227, 24.537937464248863, 14.503008746359699, 15.042418160045848, 14.92213401913991, 11.779636672274101, 12.157992709413271, 10.066971120855499]),
array([25.38116241293267, 23.29389594678201, 26.038198863802798, 24.62233074521327, 23.834903523790118, 26.276612495252028, 22.148624532138744, 20.491668881021244, 20.69731665672299, 19.85950698335683, 18.546714906450237, 19.678244885228843]),
array([25.673562656227563, 25.636297179961485, 20.587403378328627]),
array([24.49421933076009, 24.41264842437952, 26.238572860462188, 23.907211552407716, 23.70565833871741, 26.126062025434518, 25.378563010108, 23.766122393335507, 25.057521116127468, 23.64119113272509, 25.042386141071315]),
array([24.63400708032386, 23.050402280933348, 24.582544136611485, 21.17415677944691]),
array([25.749688018831034]),
array([25.526623963361757, 24.826220588329225, 25.4060231188137]),
array([24.459730010133562, 24.902808343842388]),
array([25.059722263426817]),
array([21.188013379912633, 20.940089965813726]),
array([22.093634471278516]),
array([22.84314759871234, 21.2773990031665, 21.900534825130407, 21.09624529899932, 20.628043379843657, 21.43461954510893, 22.4397144324631, 21.79419690610503]),
array([22.442749552675554, 22.5733477107605, 22.67412486182384]),
array([21.588618000628646, 21.310150417047403]),
array([20.72720382008019, 20.295474851028033]),
array([20.50510596520969, 20.60070379911618, 19.774492985852568, 18.465426527562535]),
array([20.808356086216985, 19.842419934857478, 16.153221936422103, 17.076049240042273, 19.81990566643214, 15.813418467290402]),
array([21.56538828695464, 18.492523332235457, 17.40082390596566, 20.187575265930157, 19.358583261365737, 18.500141045457635, 13.908787485030452, 14.344809154883162, 14.954431514050953, 15.81843022535832, 14.79169153248346, 12.869521907101097, 12.823898625142533, 13.262793922840041, 12.616347288255383]),
array([20.750228542880453, 18.898312687195116]),
array([20.48928356829221, 21.115051207806832, 20.468062240696867, 19.50772686597834, 18.969276798922866, 20.360182695994286, 17.74197283797268, 14.574277668066577, 14.989381055604571, 15.913094781579254, 13.996544527048336, 14.093228696866465, 15.080809552784489, 13.810381972966056, 13.332949654343956, 12.567619481901117, 12.240252740816311, 13.551854212243159]),
array([20.146819010089953, 14.287707519839906, 13.643523047167006, 11.795434043475238, 11.184343501936922, 8.17411599442131, 6.944281192474634, 4.447411813425003, 4.238187668900044, 2.76099457325701, 2.3583671239585513, 2.246047823183948, 2.242885804572625, 1.5443301635291071, 1.7925339456444533, 0.9565121178986072, 0.27908983374295837, 0.3874741307196179, 0.0]),
array([16.556816834635733]),
array([18.306070178320983, 15.74002028598219, 14.02286473262939]),
array([15.953689493045493, 15.902040849191224, 13.290315477586292, 13.152432342020747, 11.246389527783855]),
array([19.57527521715313]),
array([17.422030261325045, 16.18470016108903, 14.110404653196998]),
array([16.945581320673956, 16.381633878047232, 14.683241280175661, 13.088835208978383, 11.826155838637884, 13.438638123200095, 12.826342069221717, 8.667517460106943, 11.38491245416384, 8.862131823250142, 11.203583879978543, 5.518376620249752, 2.735868029782266, 2.8768250348900684, 3.1359963734609537, 2.4989681092077247]),
array([17.50897269991753, 15.609941363970496, 14.261342144479098]),
array([16.182106224890838, 17.366595259598093, 17.086408222551206, 15.860886207163842]),
array([17.24175896377313, 14.269410076169612, 14.380463204864947, 14.838679555334778, 15.886196017460641]),
array([16.268235250611575, 14.975915162981163, 14.462313873751791, 14.1346836413584, 14.908488004247529]),
array([17.479641144989543]),
array([16.91117592522405, 14.409654439926898, 14.957759752830077]),
array([14.530897088755042, 13.050247628072679]),
array([16.052117227382826, 17.16554922236067, 16.696917021612347, 15.731505646664354, 15.499232262606174, 12.79137052193893, 12.633891542737377]),
array([16.29732938864916, 15.207314074071368, 14.10444169755353, 15.357558204361332, 13.87325442971763, 15.821090337532658, 12.57341829182134, 12.328217059765535, 12.572643406176445, 13.028366706107475, 13.436049092049423, 12.81958809081252, 12.995327612300612, 12.360899841896023, 8.348980711608766, 11.033379231679229, 9.109377677203312, 8.976264935255985, 7.666175130435109, 9.128753088313461, 11.449504964871107, 10.156291235245178, 8.7341671958185, 8.882115802933258, 11.17482530466462, 10.423193808923608, 6.8867658031802375, 6.7826578995363995, 6.35371059820239, 7.117732588934674, 6.0860055470560095, 6.6105863392974085, 5.3714816425479555, 5.1419400988926, 4.828894059207533, 2.8842331090057876, 3.3663189776827824, 3.256356567775938, 3.351957945382903, 2.1223533943717796, 2.547530767619611, 2.1896607681038804, 0.946820326201546, 1.468375553725335, 0.8174512542173227, 1.44898269625272, 1.3059099022900185, 1.6936997987962468, 0.7533279771227998]),
array([16.45402474533817]),
array([15.70489866444513, 15.399672520904542, 11.634946434014164, 11.45883833111076, 9.51049760302385, 10.621789613093812, 7.015292490194247, 5.083421511183521, 4.602849518030135, 4.80042353330939, 2.728861692727762, 3.5350644307623824, 2.7449377023810153]),
array([12.294439615839135, 8.836669815853007, 6.98176788912336]),
array([16.08253967933058, 15.761997462507658, 15.884406127176147]),
array([15.628228408560956, 15.417756060642988]),
array([15.194192919077855]),
array([12.528265953072443, 13.662329110934689, 11.63365616337838, 7.800575948161581, 10.188501719246029]),
array([12.448924394560171, 7.5216529328158295, 9.775015410222306, 6.324612341206847, 6.353752597549527, 5.425985217497694, 1.9282495666870814, 1.8941966751745072, 2.4236462151905096, 0.6439609823998367, 0.5736936667315498, 0.0]),
array([14.182721377193815, 14.454134906196192, 14.148948530295868, 14.171299382205614, 13.480556357704486]),
array([13.209109006491765]),
array([14.0923093964092, 14.329373971689593, 11.70969660433378, 12.709285479225349, 13.70678846562117, 13.677491621169235, 13.298906494109179, 10.492255304208102, 10.893474038774036, 8.586625735866347, 2.7095831598633993, 3.449385692674652, 3.5940011511963084, 1.8338683085902843, 2.2090921212803183, 1.8027533191883092, 2.5734481559552815, 1.8020992235372388, 2.4130344333158105, 1.7329217778394865, 0.5758921261608182, 0.5579800787523542, 0.40910348196992863, 0.556528032284547, 0.46366289307518954, 0.0]),
array([14.190098956874987, 13.47962935415484, 13.78532346388597, 13.289528521110299, 13.190946044007157]),
array([13.901271055300125, 13.385340153489109, 13.327221012031824, 13.164926711663053, 13.0891538026332, 9.290491993409757, 10.781737112108175, 6.072779043879429, 5.35127328782233, 4.728172761093695, 4.895306073124594, 2.616496626058711, 3.2605144355943514, 3.5128139512241208, 2.454891293559611, 1.9164966597490687, 1.441379271702337, 1.6021849116585958, 0.8595261519604769, 1.1658594157524491, 0.7600563950654063, 0.5129028863518919, 0.5353026544393121, 0.48285287741333616, 0.0013756702437793256, 0.0]),
array([13.467965513968124, 13.15868245041126, 13.366804650183843, 13.43133836065465]),
array([13.736197790339231]),
array([13.221638337028596]),
array([12.392595651814464, 12.994346298214007]),
array([13.548472088249845, 12.14459428797933, 12.305575407435661, 13.339680316981614, 12.836710433981578, 11.65825767943488, 12.970084350825557, 11.789794609317669, 13.212869642987922, 12.657771442732152, 12.77637069396216, 13.05564490528918, 12.907196058619078, 13.596678332659982, 9.298607159659085, 10.986897672781584, 9.352574168799219, 8.166543994044881, 7.653593748072372, 9.366289149902663, 6.730323738520365, 6.383506125998786, 6.85745187021203, 6.085204405990214, 3.792438634117338, 5.197409690277352, 3.4637254023338544, 3.1382905575131104, 2.5484640551956765, 2.2884178936503, 1.892847226940292, 1.010658494946108, 1.076639009413642, 1.1942817008790376, 0.9429600732840354, 0.744310492046119, 0.2232090598428461, 0.7806055308791343, 0.13638939163150987, 0.7698746700071419, 0.7685698559600196, 0.48369231607518, 0.635273677125918, 0.12928679192184045, 0.020729955685917095, 0.0]),
array([11.952819676743855, 13.405898888250181, 12.105222031436705, 12.399633319794745, 11.742327254679624, 13.38723368820631, 7.331562878294746, 8.131839715264892, 6.8942406069907465, 5.659995502496613, 6.776373809921586, 2.758674023526434, 3.0926042677313417, 2.8434092539664038, 2.3734191017563515, 2.1259266034488387, 2.5789342289271424, 2.313992521672754, 1.9334178177400867, 2.1377946848862384, 2.231436815793042, 1.258966138574346, 1.1452247279298258, 0.7877190039131285, 0.9721896244110784, 0.964972279147781, 1.667755921447671, 1.095434948748331, 1.2176582440302233, 1.1161009473198176, 1.013239357226713, 0.3102867159844897, 0.2842264645404685, 0.5611600700367494, 0.31012534637291483, 0.5944582914558241, 0.6458114630665864, 0.2403366056656817, 0.35214193093159935, 0.7444899772322926, 0.0]),
array([13.162749381952548]),
array([12.2596252478044, 12.122342864243375, 12.074974734730434, 12.694954270251818, 12.507426774093688, 12.869421092227716, 7.321063719775754, 9.392762274282092, 6.350355451323271, 6.026378981636262, 6.93915798396948, 3.7714950749700282, 4.964846856319145, 4.71569704758674, 3.422168362168381, 2.1012616701083338, 1.8022983906228287, 1.0791640939495735, 1.6877666523588362, 1.588708432541153, 1.5172313168661065, 0.4595300188968614, 0.2827599961402153, 0.16564973784211412, 0.6377546249635715, 0.1902518366228093, 0.0]),
array([12.874287038499334]),
array([12.242866767277844, 12.039422016831917, 12.405727345577588, 11.017917018241535, 10.869950520767318, 8.677566334767135, 9.7039038729724, 7.423634331419381, 5.426477671772956, 5.724922979185416, 6.625387469337739, 4.470044926321528, 3.1716112074504133, 3.238513784803375, 2.2134866349157636, 2.3465625492011473, 2.4022281634945806, 1.809257662455334, 1.2398631941392408, 0.181385790700605, 0.3264121480089844, 0.13423948780581185, 0.646256827536651, 0.0]),
array([12.157610878143752]),
array([11.916922310308456, 6.8400681084694215, 4.939824102750562, 3.083583423005802, 3.1403046284832823, 1.8912796575221975, 1.8880188162111846, 1.166955594796622, 1.4306928431678807, 0.5393864503034237, 0.5350982217802442, 0.713614166717617, 0.0]),
array([7.912320692011306, 10.268502594284438, 11.025783366973664, 7.64264585304341]),
array([10.18098617007586]),
array([11.848644176510405, 11.830528992517138, 11.834596532073467]),
array([11.279143977660492]),
array([4.0562749721080795, 4.338393522684085]),
array([3.7287935327974107, 2.5150186970111403, 2.536225440804646, 0.0]),
array([10.549075187430281, 10.325023623962384]),
array([9.639161326428072]),
array([7.570969224808536, 6.096702506514191, 1.0406837852751312, 0.0]),
array([9.766829476229232]),
array([8.089479204087347]),
array([8.845085233469458, 9.88677124045677, 7.977806661665257, 7.076350719881476, 5.893670618216446, 6.224236293201165]),
array([9.353073345676625]),
array([7.4522162878937515, 9.132866891880834, 5.78399601778456, 6.770082242723729, 6.737386138828559]),
array([9.365459996762334]),
array([7.51188236132965, 6.048328981944634]),
array([8.300126828703876, 7.810341078923592, 8.158458080728298, 8.864271635593413, 6.83122326375854, 4.681091320842418, 4.464769596111243, 4.856790130546076, 3.4709954196815302, 3.3040846008021925, 3.3805907085756335]),
array([5.585927334224648, 3.2556545257692746, 2.28398489200389, 2.0594994458495512, 2.109223102098121, 2.0720891744512144, 1.539125927907604, 1.0743263536338925, 1.0613225412690928, 1.16637386545634, 0.6590627106481742, 0.3256987581525925, 0.17624569960332026, 0.0]),
array([7.670849941697362, 7.245135042884434, 7.032131411228561, 5.983856613117469, 4.385538433647006, 4.2349428509080145, 3.404008348084038, 2.8656592993809737, 2.0782113836463356, 1.619071739612316, 0.7657231089693376, 0.0]),
array([8.359911797908556]),
array([8.81068218570571, 7.056930049922062, 6.693080402024213]),
array([8.75656402809803]),
array([7.984714003334888]),
array([7.686835719673218, 6.1937686133801675, 6.211806446374285, 6.415935882373493, 4.131109671244811, 2.5474705044450534]),
array([7.868505977273155, 8.656636125166777, 5.442916649981354, 6.4073064831189885, 6.314915008337273, 3.0502369909521843, 3.316661831797563]),
array([5.429596943841161]),
array([3.028105870371938, 2.3728631888679406, 2.0532577284410096, 1.0331025481644498, 0.9182637257404382, 1.3658156366824625, 0.3673190103688181, 0.6493787872681465, 0.49896550060764217, 0.0]),
array([7.294570635512939, 8.368713359108684, 5.35153733624835, 6.376430668492915, 6.4038320624100695, 3.6807100779692465, 2.9001040616566724, 2.487013184020326, 2.017686024862406, 1.8325303803870017, 1.7869463734087534, 0.5117145982566256, 0.0]),
array([7.639046752231101, 2.2350605028404615, 0.0]),
array([6.19196220448955, 6.4645386224141825, 6.609765423577234, 2.989363953618644, 3.1068069912497487, 3.545457273949425, 2.7174934469029592, 2.045441172571322, 2.3339143744026956, 2.047308842004257, 2.5780583666434147, 2.047349383952073, 2.1408290382946413, 1.6078781004329339, 1.1263501317778792, 0.910983492138946, 1.3180109039062153, 0.7342946624700136, 0.29908429033649997, 0.31564043586158375, 0.37474069229612444, 0.3978511484240888, 0.4261476927087535, 0.0]),
array([7.740733327500956, 6.10823705400253, 6.455474833213816, 7.1137764582149945, 4.2772850509678335, 2.687772031715993, 3.059576350232869, 2.802412140705627, 2.2023316518271834, 1.9517298434870822, 2.056160440451783, 1.7027568992881128, 0.3457103450037262, 0.3652705957391095, 0.43583734275815095, 0.0]),
array([7.482081708477847, 7.896345120221624, 6.911856513373492, 6.865442616720494, 7.16004356719147, 4.512376150758621, 5.12380556345772]),
array([7.273069540246634, 7.089575890689288, 6.522572055036384, 4.931961691880754, 4.141421155464143, 3.915131254065635, 2.576042319614245, 2.4419164635232224, 2.520277117911347, 2.3928003542932803]),
array([7.017924667317919]),
array([7.506901677337451, 6.992342249717963, 6.676254974058274, 5.0540223112010105]),
array([4.355234189888522, 4.635615186571102, 2.9309711156637834, 2.390449560790155, 2.0684618725950688, 2.402226536779708, 2.060374844508522, 1.962452515466166, 1.9607832196917903, 1.3102085837469721, 1.0740155214667357, 1.143332993228797, 0.4706277253383652, 0.0]),
array([5.419381425862503, 4.401779605432147, 5.2290106452688585, 3.2991813514519603, 3.034112248241455, 3.101669275553349]),
array([7.546417046419306, 6.100177224183239, 6.288378534561776, 4.392724202516001, 5.32224792619102, 4.1258166653196895, 3.5155688572093085, 2.7694530689865595, 2.6556893147015703, 2.3147358375203986, 1.8625458455532242, 2.564905594832626, 2.2664295951077564, 2.143594581666026, 1.928413581055796, 2.0540657832157714, 1.8659545602666086, 0.976610482099319, 1.5244110197480039, 1.5329970072565546, 0.5795102912503673, 0.5217202425770275, 0.6149009356192163, 0.4537387113288232, 0.2615317488352239, 0.39949532819280675, 0.15771533850034603, 0.46620686040960874, 0.0]),
array([7.349928081010479, 7.962228337563293]),
array([6.9594823015087455, 2.0526017156528753, 2.2141886505402866]),
array([7.029242423733523, 5.338786640225092, 7.060621687626865]),
array([6.098440216833032, 5.119759891802081, 4.735543446996701, 2.988720298547588, 3.400228122481033, 3.2279301448984077, 1.8725594231769378, 2.4714068606197803, 0.0]),
array([5.001622720966952, 4.8995996692234565, 4.680271601771849]),
array([6.6941466188392305]),
array([7.049784835472263, 4.9974744714367425, 3.051760534497593, 2.9843574386959264, 2.253349494259439, 2.2291269118920436, 2.515716499482343, 2.2069631449147167, 1.1524808786027623, 0.7810403419880856, 1.1954417140221905, 0.17188936401789956, 0.27436489109346585, 0.586823175778488, 0.0]),
array([7.115972355377824, 4.362063929270115, 3.122708123881687, 2.1689940285123783, 2.1298204123610733, 0.14850149389706213, 0.3103838682811722, 0.0]),
array([6.571786009802414, 1.740806069142146, 0.0]),
array([6.919182278138967, 2.843477122887159, 2.0797020880868886, 1.4689665858960748, 1.6366737266170581, 0.0]),
array([6.254382688763903, 6.6707767844499966, 5.251728357933297, 5.117949200131541, 3.185240772075428, 2.788713381434576, 3.199589214759636, 1.9986218654953114, 2.278077041328545, 1.8772656411810094, 2.0282457346257408, 1.587972025140961, 1.4829699936413523, 0.9501726683679831, 0.9962755926322285, 0.40431581620792034, 0.20330163711846716, 0.5617192054834738, 0.2835172664832546, 0.0]),
array([3.1849214248626154, 0.15292681530185448, 0.0]),
array([4.999332353966581]),
array([5.451994243530214, 2.8809248286937135, 3.184109418553148, 2.9555136958726598, 1.6184159079698976, 1.4939861265791867, 0.9644519914742957, 0.6499734305755265, 0.1518721615800599, 0.0]),
array([6.257653067106386, 3.9395233102679628, 3.272087830172828, 1.83797312122343, 2.0085685406028153, 1.9856516667203854, 2.3846193426399793, 1.5792188989055849, 0.44615829975221516, 0.18493833273787097, 0.5971402651768958, 0.29657528460235283, 0.03100829466701384, 0.0]),
array([2.951834820501122, 2.772255735348286, 3.2041113627235807, 2.1283495394556886, 0.2528916587251233, 0.47764346595448653, 0.6650135188958405, 0.0]),
array([6.636514542730895, 6.376886314057764]),
array([3.4269062083142563, 3.1327717997761604, 2.4945240599958263, 2.4603769743375006]),
array([5.4816988511594555, 4.225743837435166, 2.838083677022737, 3.1524937680993994, 3.3478648869951066, 1.2407577458467447, 0.6264098857670807, 0.13124144580692654, 0.7384820284174783, 0.24074706254886924, 0.7002165786603718, 0.5595932838832021, 0.0]),
array([6.172525276675649, 4.376616761885536, 3.6617071215654784, 3.8106139493540043, 4.503405145156231, 3.1240426179180716, 2.797122707287909, 2.8489125426298343, 3.1760548580940102, 2.715896110843169, 2.450536386722651, 2.367296343973755, 1.917021878129228, 1.9915856624638364, 2.3166035535780196, 2.0041502364300223, 1.9317882305870495, 1.5358279459993591, 1.3301901784742092, 1.6156632515464948, 1.627302325321599, 1.6141520897395165, 0.25820832504233016, 0.443200784449364, 0.5155794769876301, 0.5915912468738224, 0.26123186913759044, 0.7195053291140104, 0.7277313847379477, 0.0]),
array([3.5924230216627366, 3.321332044754784, 1.7089509782234729, 0.3582359627613797, 0.4977355006695427, 0.09943995388049626, 0.0]),
array([6.057922585210781, 5.125826738617815, 3.633106723137525, 3.4212421897496332, 2.8601089643417756, 3.5985258812472045, 2.434433047203591, 2.437317767161742]),
array([3.874417666206833, 5.228012178077164, 3.679158954413091, 3.0784167821136545, 3.3498883845501193, 2.3476877236976197, 1.8575156256385876, 0.13108096354997945, 0.18475094687857785, 0.0]),
array([5.718049453137871, 3.546602232779327, 3.1250859189538724, 2.476792796011365, 2.4725050800905186, 2.313273587920432]),
array([1.9525708615243267, 1.1869401391881813, 0.6975504660887899, 0.6206357993605645, 0.36359030748990495, 0.5478536829463159, 0.0]),
array([5.309772378454394]),
array([3.403630318205353, 3.092337867524858, 3.510290154503594, 3.4814174786588805, 3.13329608140459, 3.4375161878008704, 2.2361676152968224, 2.222281364948926, 2.4521563383243175, 1.3070618213185243, 1.5740244903547043, 0.5192051873467178, 0.34441417922476597, 0.3579558550004253, 0.0]),
array([3.696817925993198, 2.60137507898055, 2.153052679917603, 2.4260869280928317, 1.819193431468038, 1.3330591345790719, 1.3351140018643197, 0.7933684060472257, 1.7008103281014104, 0.7802768195185308, 0.5060950464054543, 0.06060784063100376, 0.07211019710337993, 0.11359149667522257, 0.0]),
array([4.621766772103225, 0.0]),
array([4.384429518927603, 3.9969109710168835, 3.510323891740333, 2.7212197755865666, 2.206877464229066, 1.9113417881121162, 2.5186350988949773, 1.89064541189584, 0.9097621307762345, 1.6679887948787682, 1.207097140196981, 1.1974205783652767, 0.323335719994231, 0.5202952140079659, 0.7380912102945555, 0.7339769271267337, 0.49447664993247076, 0.0]),
array([4.571482170903875, 3.636489322180396, 5.084049258438534, 2.599630595451944, 3.4671286932477803, 2.795729661703909, 2.52140058026357, 0.49710897981884306, 0.5162528956128415, 0.4320504623847171, 0.19072401399197803, 0.0]),
array([4.365918496138225, 3.401476586350285, 1.899899894867028, 1.1143227387163832, 0.43011793253101477, 0.2388179291982382, 0.581744994546041, 0.6370487931491924, 0.0]),
array([2.7702523765423863, 2.8305993701354293, 1.6808079325585914, 1.236905038823018, 0.9767898580883494, 0.7156312988434391, 0.508250277284726, 0.0]),
array([4.822087638289007, 3.2919431448846175]),
array([3.047069727468174, 2.532471168952122]),
array([4.518711948067974, 3.5674916040329374, 2.633879619636665, 3.5251627738965876, 2.432468871826517, 2.312044115447482, 2.568642106856037, 0.0]),
array([3.3881190904715894, 2.0101232319570363]),
array([4.000717690860289, 2.685967336522684, 2.5544473702518435, 2.1899854532546548, 1.9378444011483893, 2.463752765177321, 1.540844074030774, 0.8976371937614792, 1.7148580800755682, 1.565146493078517, 0.6804311412262323, 0.48110883699343926, 0.16745753880676295, 0.7191427846661591, 0.0]),
array([1.91230379139613, 2.0315298874325878, 1.099093334018232, 0.0]),
array([3.514644207809618, 1.6571755625052629, 0.4117790601252735, 0.4048984126535603, 0.4742848861082751, 0.0]),
array([2.6978891058515426, 2.1044275624050233, 1.3372300613440475, 1.1877900317414634, 0.4188620259337537, 0.5461328877702831, 0.0]),
array([2.8166370585511764, 2.9697644231564295, 1.2983816687899106, 0.31672723178789314, 0.0]),
array([1.8175142176050851, 1.9601000107592208, 1.8419279013477268, 2.4804089777739655, 1.7484347760075536, 1.6229572682760187, 1.6594179020266355, 0.27486132701397925, 0.5181369308019542, 0.4161758397981266, 0.0]),
array([3.379160542084033, 2.617786177901059, 3.171581131883122, 2.9279148523442675, 2.960272968772265, 1.142425177022079, 1.341305398331049, 0.5881698662177451, 0.49462140083996126, 0.6943614880249157, 0.0]),
array([2.6791864383340718, 3.5811673325700215, 2.7544234421605758, 2.151497980641217, 2.3445022380580935, 2.140427319259872, 2.3180161685581875, 2.463638693581464, 2.3861112183183764]),
array([3.843912016600936, 3.4634481788690525, 3.451544368609525]),
array([3.020338226624694, 2.824417956338393, 3.4250324739622964, 2.236083251993262, 1.9526043929443113, 2.2288073441522727, 1.102550298554541, 0.3859168143549003, 0.609453218157301, 0.604418566536789, 0.07443504821792386, 0.0]),
array([3.381161547257106, 3.1872058485876518, 1.9964803558746018, 2.4437374026134844, 1.9640823644595262, 1.8314096291522217, 2.210145767584717]),
array([2.1574531519425824, 1.353721026511852, 0.0]),
array([3.710508339400255, 3.1256691859574635, 2.792485339613581, 2.439295363481574, 1.8974214367008424, 1.9416773476911215, 1.6272779227099883, 1.249833785809721, 1.3405740003012787, 1.446365781405683, 1.029427135246781, 1.4723437615468296, 0.2653814494773986, 0.2384887064193678, 0.5665035036608401, 0.014432537601463485, 0.0]),
array([3.4357236038886088, 2.2136845739548847, 1.263625566210854, 0.16416679839391624, 0.0]),
array([2.7762122985779984, 2.5698648101235118, 1.4172376249103245, 0.752121447785303, 0.6471358801026879, 0.3821530702468243, 0.0]),
array([3.271755001960847, 2.097724678130538, 2.4856559571827495, 2.0029195736814365]),
array([2.652380547104673, 2.8881587539742584, 3.025778197257885, 3.063117085520097, 2.6359298130505637, 3.41847019585525, 2.873018976705856, 3.1509027870435085, 2.5152947418231237, 2.2724702302697017, 1.8485114127863316, 2.0825769958705593, 1.8937448967671973, 2.03923093089196, 1.8746986232781684, 1.9878656954933673, 2.3395464191134794, 1.0912711191392537, 1.491448918408985, 1.0133605863649415]),
array([0.7518490636921724, 0.0]),
array([2.9429498082668433, 3.2070481188863758]),
array([3.3012109975812476, 2.6069985629862815, 2.633778570874365, 2.7349638586925544, 2.993409028724626, 2.690003949723748, 2.5023124144509454, 2.516983396618726, 2.116717739987521, 2.30055014596246, 1.8336137404352737, 2.444249166410323, 1.4997773894719635, 1.7246980224042652, 0.8965062201929079, 1.4382328761872358, 1.071503569483475, 0.3663773536008556, 0.7358682686069393, 0.7282387695910653, 0.6646284204272288, 0.7591555605626881, 0.5154257125695443, 0.29244898170453415, 0.462996121066867, 0.0]),
array([3.2785249065442423, 2.9286407135263373, 2.1427511216601665, 1.956127551057997, 1.8382766137035689, 2.08207662130586, 1.6624850307571617, 0.32051749075333896, 0.45887614132123955, 0.595025209518904, 0.0]),
array([2.9759226479671055, 2.9774893935267657]),
array([3.156682217702571, 1.059347433283232, 0.8489408407508499, 1.5217242221779226, 0.0]),
array([2.7502042198942585, 2.9730181763439036, 2.7262422165649713, 3.132110799219383, 3.2264741197936218, 3.228565639067944, 2.5301365928848214, 2.1252168958047184, 2.579748525390132, 2.079772417319246, 2.4818707744739292, 2.1127536181035, 1.5441257359655658, 1.1342024315356452, 1.7927058631488704, 1.456509192382137, 1.180192832687767, 0.6369433155206444, 0.3790567099105806, 0.3895360464100698, 0.35953851875665005, 0.04938867476749542, 0.0]),
array([3.261524574988955, 2.929426884904945, 3.1880164186366238, 2.8059530870110594, 2.108916233566665, 2.119162092584055, 1.9442392921494949, 1.8498831799674316, 2.1938551389310135, 1.2257836304336818, 1.2419147091650626, 0.7960899710061955, 0.21174730611136827, 0.42525486809267715, 0.20583927428217552, 0.7687081445052313, 0.2787324017829237, 0.7138543438846205, 0.0]),
array([0.19147690049820631, 0.27130191372269763, 0.2736116628826657, 0.0]),
array([2.159798143730053, 2.0237853495782776, 1.9737490162772708, 0.907742590012906, 1.489261651956984, 1.3068822173150823, 1.402741503977599, 0.026055191834405453, 0.0]),
array([2.57659183231408, 1.93947357926444, 0.0]),
array([2.573795068348983, 1.069216994089619, 0.9344874180543521, 0.5418205522406179, 0.5991782108482533]),
array([2.8496996421676717, 2.1221821075786957, 1.8822655684976164, 2.3043991015394787, 1.2689687346377083, 1.1687899051841608, 1.4530560123784955, 0.40204806301900464, 0.4033577595883441, 0.6744777672923609, 0.4782005366943886, 0.269526028117376, 0.5628541729799456, 0.7533160565017171, 0.7146159861305292, 0.0]),
array([2.331367109082343, 1.024269441103158, 0.9391861319779624, 0.31805517220331003, 0.2817686981937134, 0.29233982847896045, 0.0]),
array([2.784973727841427, 2.0624262826720363, 1.0319074006154425, 1.5790038952078704, 0.8082497754434778, 1.6908912402220673, 0.26690507005722675, 0.31112489664451604, 0.42670183533225486, 0.6029763961158855, 0.33867999858865944, 0.0]),
array([2.926808888638898, 1.985616665587259, 2.4870308564028125, 2.357801859741478, 0.38091411100152395, 0.40788430814718696, 0.5461220375502817, 0.5667167649962225, 0.0]),
array([2.888064235968099, 2.360114802199374, 2.438428539647925, 2.1043069800211898]),
array([2.7835509432143515, 2.8112139245916308, 0.5287187181926063, 0.3765532438542891, 0.09620120619331671, 0.0]),
array([2.1984276743688786, 2.0606864537582754, 2.5416537915980166, 2.0472867954049856, 2.514751387185394, 2.2648791987107684, 2.150064302320764, 1.6882916272613906, 1.271196785502975, 0.8536568475472782, 1.1108554348286694, 0.6137455214449765, 0.37042810896606165, 0.6201396931107888, 0.4409330451523772, 0.34684181866104213, 0.0]),
array([2.632834172720288, 1.9945200060390023, 2.2008037536974654, 1.9846191952894572, 1.8095734511211368, 2.2438519307462266, 2.33676557485254, 2.4734590248167905, 2.3446914975388187, 1.90977934330735, 2.098200178443421, 1.952578386432347, 1.2623421417381917, 1.0771311466615496, 1.546092169801112, 1.6815517690466473, 1.4250402112290896, 1.7015419702438659, 1.6058189738231157, 1.7925445192261276, 1.5078092217747412, 1.1097931353942938]),
array([2.367714394646733, 1.8760261916438372, 2.0445932244338776, 1.8405842993856287, 1.428754030591595, 1.116895317785355, 1.0532980503437819, 0.706416160585331, 0.7297151940805168, 0.5846707761572816, 0.0]),
array([1.9110800439938496, 2.2395542784033617, 1.6628523449551906, 0.481397967301345, 0.28132126413055736, 0.0]),
array([2.296123429156326, 1.2724751386311675, 0.1787865563912906, 0.6825918629923923, 0.2458117380223932, 0.0]),
array([0.5844528279442622, 0.2511587324296535, 0.0]),
array([0.494689773515734, 0.5689535081128803, 0.0]),
array([1.6010354050953874, 1.3270964818210151, 1.6855323331745211, 0.6923950467137654, 0.6311971293201883, 0.7051642755653746, 0.6745348429899295, 0.0]),
array([2.0399383478520106, 1.083215079229128, 0.8648207286596561, 0.8785032382137161, 0.7278523466788266, 0.1404130782034111, 0.37236561934503654, 0.24544158429613383, 0.45754793450929737, 0.03299604445932794, 0.0]),
array([0.22033208965633877]),
array([1.8287079414785532, 1.87786429470361, 1.4205727808768043, 0.5862173947704759, 0.39302512835857933, 0.005861216039496478, 0.0]),
array([1.7163033751909018, 0.972475453263212, 1.3062382894938904, 1.6853568118400382, 0.4592698472017541, 0.24781667162083887, 0.27907793618083976, 0.027753473857916613, 0.0]),
array([1.5382887443140076, 1.6813210987720804, 0.6101975411043643, 0.6159836671072392, 0.4309495364236703, 0.0]),
array([1.1501036979452235, 0.540731646764631, 0.7609854550565788]),
array([1.4696947160030873, 0.9815980884285267, 1.396458577428288, 0.26370248380043293, 0.2832652858409631, 0.5505176365936963, 0.1623223566306552, 0.4257447545287541, 0.5565039011604341, 0.30030803531757677, 0.4499031968486786, 0.3387127880536448, 0.6767121723440612, 0.33722538805649105, 0.32260443211923934, 0.47808267657954245, 0.0]),
array([1.1412123853984477, 1.68904100311689, 0.49145908133807775, 0.23080016647298107, 0.7745370882685723, 0.6623782438061941, 0.6586524191500671, 0.0]),
array([0.972182957396535, 0.7222631125154406, 0.3589090671990726, 0.26152067892995845, 0.0]),
array([1.0635358080062152, 1.4775385911987053, 0.8170506673897829, 0.44010162305675415, 0.0]),
array([0.9001256254383212, 0.0]),
array([0.260200508046994, 0.6894903204402828, 0.7590896987946745, 0.37984890192590004, 0.5060532724196596, 0.6337434163485463, 0.16070833056342415, 0.1597680642417454, 0.47688458229926184, 0.314358050446105, 0.0]),
array([0.8136684626695061, 0.8930306938525031, 0.45784345527581444, 0.7347725530767572, 0.5954591239460588, 0.44193162739337233, 0.19050963397649479, 0.6126537442860824, 0.6389415592653377, 0.32616201643222276, 0.27958990185439614, 0.6990499274960654, 0.1909812845922818, 0.5546919709118359, 0.14327204625533774, 0.5793312556378025, 0.07432683158631395, 0.0]),
array([0.6167623971007856, 0.44085345500442324, 0.0]),
array([0.2640067445314104, 0.4601973705648349, 0.6688868043911562, 0.0]),
array([0.652834692061603, 0.16499846533378548, 0.30313049328647734, 0.22206309350446352, 0.713049440802937, 0.17467937367703734, 0.12484812063130082, 0.02878258450857342]),
array([0.4830039268881559, 0.0]),
array([0.1932499049336256, 0.220027709348235, 0.7026247575598206, 0.0]),
array([0.4073375015007152, 0.42609415091678454, 0.0]),
array([0.5015366980197556, 0.19060185556173231, 0.13188729202669147, 0.2702535640596289, 0.0]),
array([0.5105052736446194, 0.383153314832123, 0.4100113482110639, 0.0]),
array([0.2647322040534025, 0.2687046149775483, 0.168206074174219, 0.0]),
array([0.1927694873746087, 0.07805813048977817, 0.0]),
array([0.19590514776570983]),
array([0.1968084766423312, 0.14853624651715847, 0.02364753152675189, 0.05024727140328422, 0.0])
]
d = [data_1]
names = ["45"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T13', 'T14', 'T16', 'T18', 'T19', 'T20', 'T26', 'T27', 'T28', 'T29', 'T31', 'T32', 'T34', 'T39', 'T40', 'T42', 'T44', 'T46', 'T47', 'T48', 'T49', 'T50', 'T51', 'T55', 'T56', 'T59', 'T60', 'T61', 'T63', 'T64', 'T65', 'T66', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T85', 'T88', 'T90', 'T93', 'T95', 'T96', 'T97', 'T99', 'T100', 'T101', 'T102', 'T103', 'T106', 'T111', 'T112', 'T114', 'T117', 'T118', 'T120', 'T122', 'T123', 'T124', 'T126', 'T127', 'T129', 'T132', 'T133', 'T134', 'T135', 'T137', 'T139', 'T140', 'T141', 'T142', 'T143', 'T144', 'T146', 'T148', 'T149', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T159', 'T160', 'T162', 'T163', 'T164', 'T166', 'T167', 'T169', 'T170', 'T172', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T186', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T194', 'T196', 'T197', 'T198', 'T199', 'T200', 'T201', 'T202', 'T203', 'T204', 'T205', 'T206', 'T207', 'T208', 'T209', 'T210', 'T211', 'T212', 'T213', 'T214', 'T215', 'T216', 'T217', 'T218', 'T219', 'T220', 'T221', 'T222', 'T223', 'T224', 'T225', 'T226', 'T227', 'T228', 'T229', 'T230', 'T231', 'T232', 'T233', 'T234', 'T235', 'T236', 'T237', 'T238', 'T239', 'T240', 'T241', 'T242', 'T243', 'T244', 'T245', 'T246', 'T247', 'T248', 'T249', 'T250', 'T251', 'T252', 'T253', 'T254', 'T255', 'T256', 'T258', 'T262', 'T263', 'T264', 'T266', 'T269', 'T270', 'T272', 'T273', 'T275', 'T276', 'T277', 'T278', 'T279', 'T282', 'T284']
def get_taxa_names(): return taxa_names