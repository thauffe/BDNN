#!/usr/bin/env python
from numpy import *
data_1 = [
array([30.210088020037116]),
array([19.47280950915552, 20.332260102872247, 8.652414717124103, 9.203651289088812, 6.4349931333863175, 5.782472443946339, 6.163814426793186, 5.725861181278242, 4.078964718357452, 3.074764816259706, 0.0]),
array([26.520032020617553, 15.766263093250844]),
array([25.544681288990887, 23.51944122719705, 27.099377562197585, 22.30417172689105, 22.295227255192536, 21.970136243936818]),
array([20.073592843150713, 17.935728828694145, 18.944867655194578, 19.722289015317035]),
array([25.751172510490466, 25.782103037084887, 25.570149266877106, 25.299143068168366, 24.356369415954152, 21.648823197069994, 21.580251628994702, 21.21244251370069, 21.864049025485848, 22.41616484056788, 21.24409897617413, 19.465051293948825, 17.339856024669224, 19.59706311345236, 18.565812517823627, 20.236151731375976, 13.871227902263392, 15.28899156995659, 15.89664219475879, 15.947138910066121, 14.027997831958363, 12.793207082025898, 13.166644520938112, 13.676639090627436, 13.60750856880279, 11.691283483171357, 9.704644953653949, 8.838227773719927, 10.341165004518274, 8.334574121919456, 10.163220608229771, 7.917440581591368, 8.056072046418015, 7.28724498798471, 8.97238403933521, 7.669011185238908, 9.7325079448556, 9.653546894390448, 7.841231619971275, 7.89541789029068, 10.683969044778415, 9.549574856936907, 11.608484987128113, 8.392143356371266, 7.995383251928375, 11.616129031855557, 9.371588110635734, 6.028057415176183, 6.820284805138027, 6.702899920168759, 5.641840642685891, 7.178444965529487, 5.674732784122082, 6.636658646161152, 5.427001268645684, 5.9765623496338804, 6.841671863325621, 5.176453814667798]),
array([23.203544598824557, 24.345989885620526, 21.359573643246865, 23.00301083576715, 21.90959454558969, 20.77831036075413, 22.606837928408897, 21.840985755135222, 22.73878826237855, 18.864305904818305, 19.74952006564368, 17.673472590742463, 17.089259104594166, 17.46505871421008, 16.28318568640805, 17.482547370805083, 17.364272447304142, 17.10728275740954, 17.647979962301896, 18.249197374207675, 20.358712507367038]),
array([24.133153440820735, 23.66798135386979, 24.193259670535166, 24.295840718693402, 24.221144021078832, 23.182654692638213, 21.8688223356664, 21.75817149181217, 22.47601380087894, 21.734438042052783, 21.322336044881933, 20.544968541837125, 22.88930695184763, 22.616641836234074, 21.669450985406357, 21.316398433833943, 16.90283678541158, 19.725293077399517, 18.723190343836077, 17.697179438958514, 19.980573661769306, 19.797454951474432, 19.010115763308207, 19.51723934614635, 19.025811262656635, 19.7055147662411, 16.63817007234746, 20.2929365717828, 17.573795412653997, 17.124922382530855, 19.251624866359712, 16.55993662021484, 19.554634936962092, 20.169801258694786, 14.57584871204191, 14.17822561417561, 15.314493188123318, 15.916109033778087, 14.346059228194985, 14.117618972984328, 14.913024571845536, 14.65750082537724, 15.212467028207717, 12.967136878849185, 12.621040222385785, 12.543472477589551, 13.231888725426025, 12.373254723499468, 13.54308472061919, 11.967552559420788, 12.2552768602934, 12.164144575808834, 13.602859674731306, 9.327676992512776, 11.034320098160627, 8.194328032266235, 9.266239830471351, 9.120865613155555, 7.901085524759231, 8.97143689795988, 8.859586724114497, 10.565784354864363, 10.99011908223407, 10.645403860917689, 9.405106467336543, 9.42385591556696, 10.404435142272536, 10.185195422281215, 9.936293337266555, 8.412572281993654, 10.824974269570184, 10.676868983073769, 8.936757834188068, 9.89224141499273, 8.667707101332702, 11.055892232145036, 8.18671238066505, 11.09609055799887, 10.10531251524398, 8.349380264332321, 7.148948965276199, 6.323018419048116, 6.9349985532631875, 6.737525483649503, 5.529970118633042, 6.9236869280042805, 6.506018049875564, 6.274522107183234, 6.723144073939887, 6.258715682659895, 6.267860003373269, 7.150173665710876, 5.502610529115119, 6.295107028950575, 6.3937377952017425, 6.5912428935985385, 6.7709083331853455, 5.905550121543713, 6.0298390374419455, 6.4078497651622675, 5.1067482318429, 5.096400621002182, 4.821076083303528, 4.432825512058113, 5.169215343534184, 5.241928965737399, 3.6159890596133915, 2.890400132617162, 2.738251812475789]),
array([22.6108488855809, 22.592018093962643, 22.546140402975865, 21.128923284326024, 22.821582622906924, 22.10831401832517]),
array([23.32920515373836, 21.48072494141615, 22.560655820672594, 21.54555486790392, 21.858629421608796, 21.090172431425607]),
array([21.651281147538874, 22.718879950484475, 18.10839207365713, 16.061137404055042, 19.223015196413723, 17.970638346547027, 14.847782817815133, 12.738457156735954, 12.873132905512973, 11.890979368958565, 8.014151684800563, 10.823902886107227, 8.090495379502803, 8.119805259561312, 10.841321737287979, 9.529942788315521, 9.675305074266747, 11.546667340759603, 7.814027006857614]),
array([17.732175122764062, 15.514214994264481, 13.746676851147337, 10.15536727465109, 8.936134577048081, 11.559908668187473]),
array([21.909247808779945, 21.320797728172504, 21.837186758433162, 21.447667187563013, 20.969257932977058, 19.10231221291347, 19.25007891794638, 19.23399593823097, 20.33454245467582, 17.264895405754118, 16.46918797170829, 20.108188431469017, 18.799674866355836, 15.369402120512879, 15.443647358163599, 15.300770948096137, 15.843506106610606, 15.418206897531487, 12.067982932429738, 13.02330329423697, 13.060693472431907, 13.782674634325147, 13.558390997696046, 12.630059017534492, 13.343610235691262]),
array([21.244016038395017, 21.79760720399926, 21.408849297385224, 19.615452256531686, 19.654174666281257]),
array([20.528733853838116, 21.649268310698105, 18.628762180664847, 17.411674347636815, 17.268518660044013, 18.0949503427782, 15.273370734615577, 15.799940189792606, 15.77706376978336, 15.075521709278444, 14.036326243074006, 15.565390572302395, 12.231925324514581, 10.43765282295593, 11.154933620625936, 7.89534796335048, 8.445896811641944, 9.622757927801754, 8.976299090805064, 10.330335974746832, 10.5817909234004, 9.688492384809294, 10.882062692789116, 8.601015220359688, 7.305523036217685, 6.7386875371640125, 6.917027914267849, 5.6002234613726145, 6.536273242655443, 4.783185387814883, 3.2395047487526685, 2.9877347231653295, 3.0056026074252604, 1.8348178201942904, 0.9552781324915784, 1.7798016725388803, 1.0220362323806365, 1.3407806188777518, 1.0011127131153499, 1.7690873471945856, 1.4575618079120012, 0.6117665340307007, 0.4069144996454508, 0.4810516348552509, 0.0]),
array([15.947314951697924, 15.858928737506258, 14.376214580288492, 14.278602081595766, 12.420121081206913, 12.941349334132257, 12.127615934492855, 10.963596605903291, 10.816483284399439, 10.753978383031184]),
array([20.542943622728544, 18.16917015916972, 16.456022595679165, 19.31237072872577, 18.52050076361356, 20.174855428796565, 17.472434235615513, 19.07968189443986, 16.362603438305918, 18.234558265354416, 19.179817588572245, 18.79881482905193, 17.240822131425663, 18.885499798090656, 19.040245519793324, 18.235787398032617, 20.06923192291383, 15.347272319715685, 14.886771922926265, 14.349231238178508, 13.997415145868104, 14.339582851566014, 15.442249098372573, 13.897041880520774, 14.70990198732186, 14.84505318005294, 15.829928523680136, 14.558508343623842, 13.91710081533375, 14.036989858100887, 11.881929993623322, 12.18531059884681, 12.330982585428682, 13.204340469534326, 13.723172173731474, 12.177037183337587, 11.861978610285197, 13.516152192552312, 11.247962738593165, 11.352620362717614, 11.094671110312245, 10.866815722818027, 10.841759218179702, 11.085775571070961, 11.597920980558023, 11.606073165842723, 11.457459651730796, 11.484765163434458, 11.49526902064984]),
array([19.918743850697247]),
array([16.168433599141682, 16.852858042053075, 16.690273137577496, 16.222912509978997, 18.72082742813393, 19.252945307602698, 17.505370506810795, 15.399492951819862, 13.934881594435486, 15.89793225507393, 14.392091340980155, 11.73740123692987, 12.771531233438807, 12.950072154348424, 13.816242921651247, 13.074960903554176, 13.207497672074902, 13.37918334961774, 8.95524760775382, 8.99710561560238, 9.266342042840314, 10.12501358698184, 9.26400633440898, 11.525236193268778, 7.300092366752268, 10.871996945524577, 11.276150187911476, 8.944863768125458, 7.980737133868814, 9.125591655713958, 7.940756763460737, 7.699220159534091, 9.08587018896358, 10.497017923502606, 9.569168905681632, 7.704573657284621, 8.119953745690433, 9.197142579826046, 8.744678026867557, 9.176507335705857, 6.047154486295655, 7.091634994288755, 7.024761808916458, 5.698357011484229, 6.294362520485117, 7.1723041315134575, 5.663405472292796, 5.87189045672151, 6.828997292997615, 7.135738121353905, 6.993343432196944, 6.812295935460712, 5.478240046025832, 7.208102053139223, 4.655639096980364, 5.278520746170477, 4.943647840694789, 5.128047428159783]),
array([17.505096851163714, 14.072337255484868, 12.491942174860414, 12.164695793842345, 13.697523006482523, 10.580116641525883, 11.275624006166368, 10.242722443476206]),
array([18.081337488670435, 17.933810773117386, 18.35892062036668, 18.05963619150286]),
array([17.14514017498817, 17.63312719171816, 14.631752437887592, 15.67368861836327, 15.518252700310843, 13.956926810641352, 15.36895463432041, 15.662825094453963, 14.52350692795979, 13.676602334911722, 13.506763758121977, 13.485576037070192]),
array([17.540400524005847, 18.02997201593932]),
array([16.739170314870453, 13.575447228726812, 13.433594813187241, 9.591474520155964, 7.70837271798485, 10.857595094932421, 10.13701402347052, 9.418372607556888, 9.846637014159327, 8.274981617358451, 8.531180217256274, 10.258196021425363, 10.263858046923097, 8.449770898387069, 9.267958164827588, 8.19698203111119]),
array([16.607144965850864, 17.437090646845643, 14.250013860354098, 14.616653348694271, 13.858676762260982, 14.402642734517777, 12.28725400416244, 11.680322092092332, 12.743192512936428, 12.24077894995099, 11.67550972082791, 11.679861027341058, 11.61005806853695, 11.467005271466316]),
array([17.339688221292022]),
array([15.75891149219782, 14.506542992227004, 14.024373098088711, 14.704211081164816, 9.82398982757015]),
array([14.47715904371038, 15.77612586163762, 13.146948554958625, 7.391685514647497, 10.01721259610234, 8.880194883432248, 8.509941371516577, 6.780691020288144, 6.802944241737826, 7.208838698740685, 7.201359904702969]),
array([16.37656279730249, 16.66318785066618, 15.475505573229341, 15.765821987567104, 15.961662504221964, 15.350500122706976, 15.408481061525737, 15.381562907409347]),
array([14.611280059865585, 15.092416276347384, 12.939599180314831, 12.05303840481541, 12.705791823475835, 8.477714594361522, 10.405942081860204, 11.609144621373241, 8.490038886859574, 9.233009790834824]),
array([11.757514582894109, 11.606521652238019]),
array([14.490869815552488, 15.440193330459877, 14.980945415991984, 15.400817422775411, 14.217141306571715, 14.593050839861617]),
array([14.60248067932399, 15.348745226758403]),
array([14.580460468027464, 14.150586882436421, 15.032292372038254]),
array([14.636431968584482, 9.362666029790283, 10.8502947642339]),
array([13.88854370105447, 13.889611599974799, 13.477169769344746]),
array([11.698423035000573]),
array([12.112852851402277, 8.256484728855325, 11.327511964978074]),
array([13.422043581007397, 12.837803567124256, 13.522776831463878, 13.77736886935351, 13.645985408864664, 13.360080957650835, 12.918491392725565]),
array([13.737218675073018, 13.117986848723968, 13.153913887672209, 11.308121487239681, 11.396596363833284, 11.400915871497356]),
array([12.52671282083167, 11.95548478770878, 12.965989325554734, 12.055845301035774, 12.048241821898003, 10.344901591226524, 11.451682372703514, 10.933260264038225]),
array([12.227730833101685, 12.76901118078148, 11.64451460505484, 11.045886817004783, 9.346668326632381, 9.35969533781483, 7.676268797546477, 10.963502939783936, 7.998664484232556, 10.63387305636222, 10.648463978395103, 7.201164114129591, 6.811948762696107]),
array([12.729461407618514, 12.07227784925765, 12.74731277781721, 11.685848039228885, 12.370491525471763, 8.603781373941455, 8.807123512784532, 8.999123575016956, 10.180915906885682, 9.779354355617743, 9.180165151659484, 10.890105757882683, 10.178852240037296, 9.508014740110951, 11.297434870642586, 8.302969831564479, 10.511945047425803, 9.209505566235677, 10.731429533864757, 7.506609855459388, 11.065056398027668, 9.925509553254129, 8.166474063806398, 10.27039339144768, 6.2814565186158795, 6.4592261411340415, 6.948528386595238, 6.225190269668242, 7.120992669002647, 6.599643431351517, 6.6297965235133764, 5.622205147068759, 6.386788588753478, 6.480288166693998, 7.052638449638745, 5.568602246170606, 5.918349396599717, 6.6527378598501965, 5.776748478632896, 7.203537095070844]),
array([11.77540465197785]),
array([11.649797762320174, 12.407876592204133, 12.212245768326705, 11.694063476610877, 12.008665650634635, 12.355978281947577, 9.736080626555607, 11.278510608256578, 9.848927130591717, 10.21585464510409, 10.188090633833168, 11.00812620300383, 10.860411773428453, 11.430382949177353, 10.780720342101626, 9.786892833132612]),
array([12.12238916746217, 11.539301118706614, 11.27600006326947, 11.07087829552336, 7.368419218321842, 8.03038792237951, 7.567042356086212, 9.752411497421347, 8.794405562657257, 10.630428731654575, 5.80637188031618, 4.8380676356767545, 3.9614666581788596, 2.915576547059385, 1.0144143682687667, 0.9680463190584376, 0.0]),
array([10.014781661188996, 7.7777980031661205, 4.254561471664853, 0.0]),
array([11.313143973542948, 9.06966721777927]),
array([10.662297803262325, 9.27326209200752, 7.910123033180288, 8.394765964357944, 10.863849717803587, 7.536050097949964, 11.021346622520568, 9.534151413499611, 8.10698724245827, 8.084441912993967, 8.405659728590084, 8.576553154949053, 5.553496499443121, 5.608235272526063, 6.549487323469004, 1.2562398911358696, 0.1450696854559016, 0.0]),
array([8.89001763006864, 7.453221973429233]),
array([9.603585408004026, 7.452615157849285, 7.380993068661517, 7.485726816180768, 7.539644406707452, 8.719637297293092, 9.618333862581638, 9.80870718306759, 10.927862691461932, 9.916996276278201, 10.254879265370437, 7.331586891113124, 11.028251484895478, 9.526986193524086, 9.49123959639878, 7.583926606921637, 5.4382813581715155, 5.8135596541790875, 6.272770421044554, 5.740001759861905, 6.566950757867064, 5.564105665398067, 6.555383893451575, 5.5132683374544325, 7.216741640276689, 6.062679493681026, 6.412396640195869, 6.059450384649349, 5.769831694602048, 7.139961291121181, 6.903392080945729, 4.160727343031494, 3.73380056343667, 5.149673820011895, 4.087154190223883, 4.44049889854619, 3.334909489658558]),
array([8.76560481245192, 7.751886585852244, 8.61862330783177, 11.307503156980419, 7.502896318445454, 10.389729535437812, 10.939217086564753, 8.783483454824987, 9.386848656636381, 9.050825387588157, 9.00333182747067, 7.5142895994604695, 9.617791079408748, 9.281330377902263, 11.237079903536417, 10.498107077686251, 7.540038457201016, 10.156635923616781, 10.480028268453921, 8.784138422252761, 11.039152320510922, 9.768542074528163, 9.205181631257934, 7.5091949383431, 7.899851637762878, 9.718852715676878, 11.133205027261042, 7.673425614422407, 8.11876814845359, 8.633961463738785, 7.95716973020704, 9.044549773986345, 11.365601697931025, 6.820934228879985, 6.218868199477026, 6.240164901248251, 7.219877599306882, 5.4105859537254135, 7.176065824722526, 5.731500494676284, 6.208286686923123, 7.127705166633691, 5.473669333979194, 6.7650500828841755, 5.887270328651061, 5.953600093619734, 6.227958331718311, 4.482027670400798, 4.684971790182807, 4.538150488452646, 5.302774026519013, 3.8906128764062284, 3.947886631613188, 3.8389740949867255, 4.998612285329335, 5.232054758722701, 4.638504725229623, 3.4933016336696543, 3.288548859529646, 2.9028417648890685, 3.068699628621208, 2.0541715738259207, 0.8772648886624721, 1.2019725475574328, 1.2806004340055597, 1.796633928507644, 1.4296891812000259, 1.232596979344307, 1.1179580380960243, 1.4252432461485047, 1.1546302308278864, 1.3039396328172963, 1.7136308329960517, 0.5230194568337414, 0.7377331507252739, 0.5554047829105296, 0.0]),
array([9.358966122085395, 10.328520355976309, 11.193482494505096, 9.629387101919175, 9.204670554212337, 8.729780765950398, 10.012181738160564, 10.618177326411118, 8.933420592865412, 8.266222599076432, 9.273368175014083, 11.342589010913937, 10.068218489494956, 10.238698812482866, 10.25114698888855, 8.95333022574092, 8.555453614344453, 8.26736744343226, 10.33758606254914, 10.906280245649775, 9.188723701781752, 9.009029220253744, 9.710779346266484, 8.3829003487289, 8.844003373467796, 8.190205480754248, 10.292202957936093, 8.865230660498467, 10.274042774241083, 8.764683225637402, 9.045952450972852, 9.73298696161456, 11.143705205254983, 9.334469372382697, 10.843004498268014, 8.293989165975537, 10.812031082868167, 11.352763770718314, 9.20574741792402, 9.035312232325818, 11.531091263975211, 9.705269648201416, 9.423238116773666]),
array([11.276117653624828, 11.182982950430016, 11.059654804345742, 11.1200508659151]),
array([10.322076983528135, 10.718987715754116, 10.54443058991626, 10.165421807892224, 11.199228480333426, 10.812046629707847, 9.907829201302192, 9.900783129372094, 10.097279290457955, 11.095866096176197, 10.55076000430555, 11.083090247115049, 10.897087628858756, 10.930573847089157]),
array([11.178536536839841, 9.924870117015239, 8.280292893162125, 9.709102559497342, 7.249383262914975, 7.974942485479314, 10.796523829781817, 9.785009012957547, 5.843985665740404, 5.5362942874352985, 6.262135878467162, 4.73914230265672, 4.797595100800622, 4.65753606458008, 2.9225526362642738, 3.5279505075699396, 0.0]),
array([8.074538524249053, 5.604789272519228, 4.371434332549198, 0.9839288465529012, 0.0]),
array([10.94586435448127, 10.623938206827088, 10.54913025547711, 10.518162353070757, 10.201207123599616, 10.04124184117481]),
array([9.107977231235052, 10.64928069124965, 9.357612823634293]),
array([8.20787071489785, 8.183854352194293, 6.344287998479204, 6.593604622010512, 5.963138383095735, 6.119501718374249, 5.857396510950739, 6.682955828397075, 7.057384890673649, 6.929922982503637, 6.958420474993954, 3.519073759264303, 2.407212125199682, 2.1681368736247566, 1.34226461985647, 1.4734906879232308, 1.7900505646703488, 0.0]),
array([8.495936176506817, 8.94033999637784, 9.619608164499695, 7.391898934600128, 9.603926921677825, 9.539709322082391, 10.33905883512886, 7.315317455520148, 7.971875792864466, 10.076133475166737, 8.914961790478472, 7.7341654682489915, 8.83708267463454, 10.608418710530234, 8.680510617096889, 8.908832172123685, 8.223229844900175, 9.223951144828241, 7.452421516254083, 7.888453879258744, 10.222450801917416, 7.9401948128161655, 10.377716101664149, 5.702086912636028, 5.65933515393532, 5.412329461898582, 6.511568760191322, 5.43193380081126, 6.59021951630414, 5.833624129311377, 5.875006820892704, 6.24156741411381, 6.708997646528623, 7.0005126053440145, 6.060671793114583, 6.453885597094882, 6.299261740740929]),
array([6.578799592671753, 6.537060998871548]),
array([10.236419773398977]),
array([9.141526803114601, 8.713685357112334, 9.097799602877886, 9.112550102651054, 9.656876744948546, 8.632512476895162, 8.731504561834013, 9.27502953921204, 8.837362697231624, 8.677880329179988, 8.83586332736772]),
array([8.384945299727237, 8.94385671633755, 8.284342224696967, 8.33651662410969, 9.08766114313949, 8.865970876112902]),
array([7.9429102896103405]),
array([8.471618346930034, 8.66863553208643, 8.823188560987091]),
array([8.22278273120704, 9.028951324587888, 6.175111192282716, 6.322610161489603, 7.138964519367621, 5.213260623209107, 4.411600000811835, 5.200420090373538, 4.4784364027611065, 3.5292843976655646, 1.7145000002020674, 1.1900103531082955]),
array([8.246106514485552]),
array([6.518530346776895, 5.882645698905792, 4.497801928364488, 4.6565522093497735]),
array([8.495711882683343, 6.835111255995573]),
array([7.7319108622337485, 8.16412997956654, 7.950131822223394, 7.735297160937254]),
array([8.135107049926853, 7.383632062287528, 7.958031325847079, 8.101131528963512, 7.93516937087864, 8.11377480925333, 7.3520184135694, 7.281373532122817, 8.14444390751545, 7.334655697585906, 7.279353343239784, 7.956497467622253, 6.927022775083961, 7.002006909488048, 7.097909483903687, 6.744247132941769, 6.6492088788570936, 7.136667401320924]),
array([7.953341197515804, 7.880021329700246, 7.615312889290573, 7.708114289621544, 7.969235060374443, 7.394514472348188, 7.687806759412985, 7.991952908367542, 7.703443396960651, 6.17977416202296, 6.610460403796968, 7.220321098496432, 6.975732550400564, 6.3937701565946305, 6.158714193155471, 6.931112588333487, 6.382704110780984, 6.043421794775463, 6.5993026446148715, 6.564752938110752, 6.652841572847075, 5.896512886324643, 7.022299474068728, 7.1307598459037775, 5.979276431329994, 6.914187474777215, 6.056764811030182, 6.6442325389948955, 6.469307168261852, 6.589851023696214]),
array([8.12315914037965, 6.903033256087932, 5.43757200442366, 2.2032457059476385, 2.451445975689448]),
array([8.166083278577467, 6.169132154491116, 6.212568524467851, 6.4707757545559, 6.630025761661678, 5.714332579719075, 7.064304890216869, 6.843546559381843, 4.561938808745005, 1.8126467759871436, 1.9961700385748242, 1.650195872537088, 1.6472207064679798, 1.0734640676595932, 1.2273387938638494]),
array([8.013254601610612, 7.9929865424618205, 7.855259900928296]),
array([7.8736635889017, 7.6775867950146965, 6.2484473418862, 7.037884492326214, 6.749018454399919, 5.933180856313564, 5.706122414087356, 5.797576436899353, 4.703798302760158]),
array([7.1121914842985445, 6.9090284254457766]),
array([7.443096934748402, 7.704850517535755, 7.817382119254452, 7.445708044203911, 7.275603500726839, 5.531109827593216, 6.268588180037152, 6.766354579278766, 7.008029594906987, 5.433539157918136, 5.644658252384382, 6.26428188430285, 6.4880414651195, 6.698314027429864, 7.088253374447765, 7.019335551262849, 6.165361491652773, 6.944863256475704, 6.661860220380618, 6.078674182341647, 7.188886373075931]),
array([7.929532414287976, 7.850121636887074, 7.253604755492755, 7.64856207122164, 7.74575520207076, 5.642756909255491, 5.54377094865154, 5.906103112270459, 5.660955041504344, 5.635010933481913, 6.663809921228058, 6.531234691644704, 6.699702978661168, 6.186274862786401, 6.651999306674915, 6.855602839544422, 7.060547794550893, 5.921895596426053, 6.252810926276089, 5.550315820152824, 5.943155562952708, 7.031685030072868, 5.128086163201038, 5.231927360120591, 4.8889366340732865, 5.228692060760156]),
array([7.702470171452283, 7.5712983169522765, 7.183563836291038]),
array([7.4278322303988435, 7.215217541289014, 6.622376573500359]),
array([7.261791226877677, 7.413241420104399, 7.4892119097173495, 7.599219315174361, 7.607710466918713, 7.514531932753891, 7.27161898461042, 7.353994075323542, 7.460268048206814, 6.8738305678413365, 6.878172139057934, 6.89567402016718, 6.862848495207367, 7.03270992946288, 6.886313849588751, 6.911297154554239, 7.021956158597175, 7.192093065530744, 6.910834418328097, 7.147410949842161, 6.882838973987856, 7.191660509430819, 6.928258250210028]),
array([5.887054683598399, 6.725553224376318, 6.470979912773426, 2.154155538285271, 1.0951936946519045, 0.9867695541209927, 0.0]),
array([7.233697694458657, 6.106559003869914]),
array([6.976963938587336, 6.772082019415491, 7.176545572310917]),
array([7.172380409872014, 6.482305731852503, 7.075062657426606, 7.063461088567318, 6.081879912670513, 6.814085720862605, 7.032749776745967, 6.479537899289183, 7.068188196761169, 6.881179087437627, 6.615562031891809, 6.074970927593043, 6.800924679069861, 7.075552923198472, 6.666275219851445, 6.022324322625288, 6.525183112145699, 7.209855939873342, 7.000805842686544, 7.09161786184279, 7.016439012447549, 6.491638250551129]),
array([7.279215223088778, 7.264164900737277, 6.101778977895877, 5.719805635572987, 6.028975583809892, 6.273289567175749, 6.669050773804848, 6.881275200756662, 6.500145836632609, 6.698179387485559, 6.895325857349031, 6.026268597665753, 5.801462622122769, 6.213601391139823, 7.13063899432407]),
array([4.249032455428285]),
array([5.711822120927437, 3.9056740117188546, 2.592235072887831, 2.5549429142624613, 0.0]),
array([6.25688506835424, 6.583764846615032, 6.370408865386215]),
array([5.86814326546364, 5.67425694826698, 5.640456966556099, 6.211065261034261]),
array([6.5757051802713296, 6.482797444206608, 5.506400339075176]),
array([5.508087108538017, 6.106937843095548, 6.3193652026591325, 6.388412130175089, 5.95906267499997, 5.398911462112764, 5.440052717089471, 6.333667490151626, 4.543756215114236, 4.609631695603073, 3.708215239623396, 4.013714725612127, 1.4965053163898223]),
array([5.4073040771203225, 5.83081041522604, 6.102150168517007, 5.551951928834526, 5.574876117619926, 5.408667771664509, 6.017992691038069, 5.020350563750744]),
array([6.049313661122832, 5.9524559576476, 5.50919045932138, 4.775787778649329]),
array([6.040363058175245, 5.775078342489219, 5.4112560010857145, 2.986665052105935]),
array([5.893049165991462, 5.460550981429842, 5.783058499808731]),
array([5.555564747807184, 5.202346253907473, 4.792672467711017, 3.5541997478749945, 1.682149722319774, 1.2794502468779965, 0.0]),
array([4.549328403963656]),
array([5.525487369392486, 5.38776293374239, 4.9300096244462335, 3.6915003916136784, 4.36082468644282, 3.44889847689374, 3.4376046147128725, 2.6170093788378344, 2.5737411153256904, 1.370937999171392, 0.914358562915551, 1.425294531992598, 1.7667250614645422, 1.4255519519806739, 1.570642973531288, 1.6177205839125248, 0.21526518924127602, 0.7427914590857385, 0.5112194397831986, 0.4548029243672564, 0.0]),
array([5.435406829341249, 5.556733084143541, 5.370115826816735, 5.530605092634836, 4.10317403172421, 4.057247328712587, 4.858865901952565, 4.59728485023009, 3.896375989349317, 3.5426602109888097, 3.5558116802966486, 3.2169852198213014, 2.8725678142739595, 3.0833937102657454]),
array([5.67425240046411, 5.546441374368448, 5.714221089973527, 5.780191771865093, 5.678699375267471, 5.625182591012438, 5.424965877150502, 5.75941624512694, 5.807790662715363, 5.524528771701771, 5.5236219880916915, 5.459606442204708, 5.6568750461786825, 5.71758533804725]),
array([5.388193018983213, 4.967793093018881, 4.773868113588442]),
array([5.560840802332071, 5.343770328642421, 5.684221132735954, 5.526643954749342, 4.9874220824694495, 5.263449551775097]),
array([5.096355904787173, 4.443196360782454, 3.4467949079964035, 3.091189247675003, 2.8765535172251675, 1.9186351067647378, 0.8312826006567228, 1.7316508649583096, 0.060897357289731294, 0.1088628733038069, 0.0]),
array([5.412487420256383, 4.922157417171083, 4.963758156490055, 5.17423569168979, 5.2415746909642085]),
array([4.7283812037736475]),
array([3.038686947901271, 2.9001337010745223]),
array([3.3685220277694397, 3.1670249581472047, 2.754763525425302, 0.9003498875611875, 0.8162525813094555, 0.0]),
array([4.587536803292591, 4.002777334126554, 4.100043158096655, 4.858885861835681, 4.906732001868328, 4.724529726432587, 3.7471440980231385, 4.505401411997611, 5.037187349369747, 3.7727152591378665]),
array([4.8503655849688165, 3.7331012948021893, 3.8205738997889664, 3.7961690112615782, 2.839585649553721, 2.2050035220987305, 2.4294608494959724, 1.141362290405811, 1.0450748525386249, 1.553318876306859]),
array([4.093158225539574, 4.25440700562332, 3.649064936330766, 4.266438648142379, 4.689439752347284, 3.2329710726901757, 3.3105873164178976, 3.5526046823272677, 2.0909443698909937]),
array([4.8737533257868995, 4.721510788338054, 4.678704147468973, 4.872667186643349]),
array([4.871234039780418]),
array([4.721384338080343, 3.3206800469567934, 3.0487761035475245, 1.5515812609928499, 1.6585681360576463, 1.3366673943623175, 1.373528808057934, 0.8746015623544282, 1.510401781058695, 1.3435313543770486, 0.8151001825806266, 1.7141084972725387, 0.6171905831929035, 0.49952492000390575, 0.025938760200404504, 0.0]),
array([4.034915434541948]),
array([4.56618808475976]),
array([4.457975807330848, 4.0945732059895095, 4.494528315922643, 4.46790355433461, 3.946127865332987]),
array([3.5163048615381847, 2.7029543842318566, 1.5771640486222034, 1.588063199709817]),
array([3.6611029023258537]),
array([0.9390614196557657, 1.6258785314716322]),
array([4.16908003772111]),
array([2.1261680381733514, 1.751742104831923, 0.5976264332570637, 0.0]),
array([2.7758603488085773]),
array([3.7338326305802734, 2.9578769436835435, 2.8040732861487814, 3.3957836678156483, 3.3360175097276206, 2.757004038838143]),
array([3.7676571992433296, 3.3573035229427326]),
array([2.9338483711318415, 3.0161389634405102, 1.9489595991964448, 1.8335870013676994, 1.3647786338250307, 0.9473877255102167, 1.0524169190865575, 1.4088677935534464, 1.2633543032954742, 1.5780801090459686, 1.7717668535067905, 1.3522805632931298, 0.9285642273206625, 1.5924256270942612, 1.5411101622354515]),
array([2.674904387307938, 2.434467337259463]),
array([2.812418522963984, 2.4884790250783135, 1.005018679142972, 1.6180052884685674, 1.7238641524209177, 1.6140501571262869, 0.9119209936721561, 1.0079025436809772, 1.739977869599988, 0.8925422979521801, 0.8505232998749253, 0.4956231628367721, 0.7692681814987313, 0.7273032957980012, 0.0803370423000696, 0.0]),
array([3.0893643486816122]),
array([2.939416424579826, 2.6619329676947983, 2.700346817675233, 1.15944586471725, 1.420997217990513, 1.679333584193197, 1.7216880107451387, 1.2670899347745164, 1.1360107226739364, 0.8617685202575118, 0.9967087419337629, 1.6417588789830286, 0.9549585615516687, 0.9917813135606323, 1.5199926632332754, 1.6732233568992128, 1.2607640667481443, 0.9114980467471463, 0.7634313787832119]),
array([0.43279609619680287, 0.0]),
array([0.008308453420098422, 0.011743158554629693, 0.0]),
array([2.7472996450801768]),
array([2.6510358748109057, 2.5842252122698723, 2.8592353081827495, 1.9739481503729883, 1.1927212077455682, 0.890669096042098, 1.7564568462983836, 1.5254799835844337, 1.4774126273499197, 1.2863722906602644, 1.38332841021083, 1.0737913790962956, 1.218023824648593, 1.6736851420785084, 0.989589510984939, 1.4160004312921406, 1.6771374255786635, 1.633355802864954, 0.8716942996378085, 1.561907998989191, 0.7362721830604417, 0.4828046442460034, 0.2915108151473926, 0.0]),
array([2.688639174337903, 2.0519954392870194, 1.8268387996228113, 1.174645701598878, 1.7184991643879057, 1.6562069140294835, 1.0054682943934616, 1.6201612158235397, 1.0976932376447102, 1.361200073832592, 1.6584684019125773, 1.7776347753192554, 1.0109167098225167, 0.5054463983988251, 0.45304317992765725, 0.05926159369610552, 0.11433114420530789, 0.0]),
array([1.6204995797217265]),
array([1.5003991813123025, 1.27370349642658, 1.5688287634093934, 1.1354494575457281, 1.4912771598762615, 0.5841255324435812, 0.5014465216669817, 0.21492449679323056, 0.0939915567190468, 0.0]),
array([1.0927757379172451, 0.8665208955632728, 0.0]),
array([1.1776460178544708, 1.0105515812321832, 1.7240621874103421, 1.6811732057955338, 0.6845379581343027, 0.49020746358648914, 0.08634410481128152, 0.0]),
array([1.632462490326135, 0.532216426297274, 0.0]),
array([1.0715296160469867, 1.339205337855951, 1.0027337697299714, 1.1352114312338475, 1.391488579672632, 1.7467964480910363, 1.4256632310731343, 1.389452956646695, 1.6461507033980622, 1.1968027861065411, 1.4221827462224854, 0.9372151134745693, 0.9918317733798592, 0.6192606215638743]),
array([1.0828829785169365, 0.9558297849308961, 1.3979462454467875, 1.084396199824988, 0.7829428700946299, 1.7601002629804827, 0.43207505505170746, 0.7105632511848695, 0.3120492319059642, 0.0]),
array([1.7845697096091286, 1.587432269692532, 1.6973556028859462, 1.4238388923704786, 1.153164430027835, 1.468386202240677, 1.7759311534314384]),
array([1.7987767278348288, 0.0]),
array([0.967208536831149, 1.2983426296443725, 0.6520708439239523, 0.0]),
array([1.0140003747684097, 1.2491115342711203, 0.13628663260275864, 0.010747910531962224, 0.0]),
array([1.1301998607305017, 1.0600450723659198]),
array([0.8366402843303488]),
array([1.3602323392449718, 1.4105190579132587, 1.358640374564019, 1.0775888659992485, 1.2848832218405686, 1.3980623203200988, 1.1020190278504023]),
array([0.805609780115225, 0.27925393622137995, 0.4576911866661661, 0.0]),
array([0.7857438582977542, 0.0]),
array([0.30625900729950406, 0.10184737602916649, 0.029573439367718363, 0.0]),
array([0.15975376730177304, 0.0]),
array([0.13612225273834216, 0.0]),
array([0.04719391000802024, 0.0]),
array([0.10635479554227044, 0.0])
]
d = [data_1]
names = ["57"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T2', 'T4', 'T5', 'T7', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T28', 'T29', 'T31', 'T33', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T47', 'T48', 'T50', 'T51', 'T52', 'T54', 'T56', 'T57', 'T59', 'T60', 'T62', 'T63', 'T64', 'T65', 'T66', 'T68', 'T69', 'T70', 'T71', 'T73', 'T77', 'T78', 'T79', 'T80', 'T81', 'T85', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T93', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T110', 'T111', 'T112', 'T113', 'T114', 'T116', 'T118', 'T119', 'T120', 'T122', 'T123', 'T126', 'T127', 'T128', 'T129', 'T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T137', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T149', 'T150', 'T151', 'T153', 'T154', 'T155', 'T156', 'T158', 'T159', 'T160', 'T161', 'T165', 'T166', 'T167', 'T169', 'T171', 'T172', 'T173', 'T174', 'T175', 'T177', 'T178', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'T188', 'T189', 'T190', 'T192', 'T193', 'T198', 'T199', 'T202', 'T204', 'T205', 'T207', 'T208', 'T210']
def get_taxa_names(): return taxa_names