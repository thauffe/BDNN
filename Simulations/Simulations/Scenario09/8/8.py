#!/usr/bin/env python
from numpy import *
data_1 = [
array([28.728213955191876, 31.627280818853492, 28.780237267812204, 31.42221767317122, 30.936897163429297, 33.06913976076759, 32.53044979298503, 31.146881874345308, 31.129641386260904, 31.637980244160705, 28.252262541928662, 32.50223608157436, 27.20749976922646, 25.84850135556772, 24.08423902108808, 25.741811697530018, 24.54291616432363, 24.766950230305135, 27.68103734110987, 26.1330272300747, 25.427220873517356, 26.764965803369275, 22.296731339033283, 22.749501903409765, 22.302450069791888]),
array([30.159755073234482, 31.253616394025393, 31.278141394387248, 29.922364747425036, 31.572433714790453, 30.266876338112915, 30.706420102955992, 31.76835743351719]),
array([31.712217556219006, 30.624370551352964]),
array([28.518624518595438, 28.429068022259614, 29.888264526711392, 26.813566902174585, 26.758427555220788, 26.517489590854584, 26.655608814053064, 26.63983361056143, 26.024909688198502, 27.910955843248235, 26.31310715529585]),
array([30.459733103621183, 30.163306019667154, 29.749891113394355, 29.90096839346941, 29.783234108724233, 28.725652206835022, 29.869203184772427, 30.563850476745685]),
array([30.026713864340604]),
array([28.988352063480136, 28.60723913717658, 27.32982362570926, 26.184760170867534, 27.08558972785912, 27.229629165259286, 26.51269564899483, 27.2840381381071]),
array([29.241645445635076]),
array([28.99805040880301, 28.865290469021847]),
array([28.31851898575436, 28.4757971448319, 28.731450646346662, 28.860554998808684, 28.513920445888534, 28.11077957100658, 28.887706798758227, 28.50960781049822, 29.53237721193446, 29.423055909769477, 28.855256879508698, 28.116753361576595, 27.647291419594858]),
array([28.45801854867475, 28.985912040331463, 28.994799959805498, 29.14420104832595, 28.247659764150004, 25.471832839559237, 26.861100231340135, 28.009977420273895, 28.045223156632833, 27.568229423279167]),
array([28.47759879147204]),
array([28.11309021337558, 27.99827895001206, 28.03173954549263]),
array([24.4314474867828, 20.607499016502285, 22.756329587868144, 22.894308879593126, 19.481288023734606, 16.132414807290502, 20.34010701078613, 17.95157112472213, 10.124500306249486, 10.552347276852252, 5.518629459120137, 6.267062052306532, 2.990754007377157, 3.2390510182564745, 2.0909668934601893, 1.0587588193896988, 0.0]),
array([26.58875619833902, 26.789883087887585, 26.592982666077308, 26.51609862907131, 26.97624577176443, 26.471594259582968]),
array([24.14531448749564, 24.543913015664852, 23.536898568499737, 23.27313942981844, 24.66470367901551, 25.82721468366304, 23.17385047596717, 26.20286990527961, 24.70894772063899, 25.51005239672891, 24.472561286222568, 26.563952802718298, 23.976183517974935, 26.34556180125293, 26.136806716440205, 25.76348862303168, 26.735600898235575, 24.144444643044892, 24.168339137790163, 24.872282979011473]),
array([26.453479678167753, 26.24254591133457, 26.512847060696156]),
array([23.245452312626522, 24.282283445071837, 25.084834012506228, 25.75801039506085, 24.619183440143498, 25.770664187867, 23.31855971144958, 24.86659170777542, 23.390727958256306, 24.86650933641859, 24.250092855958474, 24.683866137533336, 23.257324709650234, 24.946570196195033, 24.42503911337719, 25.86808424198297, 24.768342161458, 24.569740154727345, 24.370820588388103, 23.577398077613065, 23.52578304401925, 24.959033325524615, 23.109872913911545, 23.210720534567884, 25.37689854342638, 25.62805449366169, 23.110943360698784, 23.041199367950014, 25.742308541735028, 23.31763755744292, 23.98586594101297, 23.039000504800633, 23.56247876558671, 24.825677380892955, 22.151028373054267, 22.92562189327326, 22.937548440361226, 22.560538492434265, 21.7596271984691, 22.440921036257897, 21.677745568327115, 22.881307926407747, 21.65505399282169, 21.929429668582184]),
array([23.959184042431588, 24.997751451678166]),
array([23.088526960222115, 24.23026415224301, 24.25215286083731, 23.614141816341782, 24.507772570588735, 22.700590847422372, 21.17022819519362, 20.68071403705796, 20.829889151730622, 21.521499305186364, 22.294990247126663, 21.83960553639129, 20.62468985978452, 21.879025977082613, 19.041233963532605, 17.698669383711767]),
array([24.20815230066185, 24.114555333514133]),
array([23.24735194564979, 23.081118083385118, 21.837550895469477, 22.575636018101903, 17.17994469439363, 18.94836705656527]),
array([22.355387695834995]),
array([20.82080567905923, 22.26905991416539, 21.569084354979598, 20.661183021828354, 20.824839560074828, 21.63037009140383, 21.716018246224085, 21.166241216157573, 21.210437953783106, 20.613492650296905, 21.331348754396597, 21.965360481673812, 21.645359533831893, 21.05771950333009, 19.22507543650073, 17.39962135007528, 19.115047141487736, 18.41771686914686, 16.905645922761543, 19.245107420021185, 19.93075084213979, 17.78759035656378, 19.43183745638078, 17.366600031982713, 18.265002722605477, 17.348813869052613, 20.09820347881651, 17.978125560941105, 16.844456436016877, 17.12172991756183, 18.900182748236013, 14.975095940849469]),
array([21.447610570582416, 21.689363428246786, 21.176304934412663, 21.30731772524337, 21.26196657341159, 21.80319389053465, 21.30798691596884, 21.326653565346003]),
array([21.644596208792148, 20.890108501348028, 21.554553478572767, 21.58296689825592, 20.539877243929958, 21.350962576696805, 21.428271016346407, 21.162114819349373, 20.87102489925501, 21.216714995196, 19.269174791199607, 18.937886431322, 19.772544576139424, 20.30163530848884]),
array([21.399510022607988, 20.70664744308412, 20.63090035810878, 21.3496045315528, 19.82739121605978, 19.854500271507334, 18.729301288488287, 18.692553379622492, 19.032407044485733, 19.053291644784554, 19.9294113019989]),
array([21.35394764244885, 21.521682153868316, 21.35025395274009, 21.36867561667699]),
array([21.261618304731506, 21.070893935520637, 20.662490736889698, 19.495815005234494, 19.145998017753804, 18.97215875519169]),
array([20.566684825240525, 21.365420039910017, 20.9925917276238, 21.127874743031107, 19.82519464780558, 20.316253717659883, 16.80482451628252, 19.7553206424469, 17.05942200205078, 17.080496897904922, 16.419823721610165, 15.546793422579443, 15.114805857297222, 15.495285115105894, 14.859200873694308, 14.23528348393658, 14.502551724688884, 14.637961872993825, 10.847654641631953, 10.672266191611401, 9.07095085591631, 7.290717697909409, 10.471347132111882, 7.696842237969411, 11.450651221724227, 8.715380860735438, 11.482195769969925]),
array([20.8914652196679, 21.00513493608825, 18.92293350574515, 18.844122331968432, 18.998399203487196, 19.024658709068536, 17.136720935852797]),
array([20.64391898503819, 21.104518006202134, 20.797546441029873, 20.769402521404125, 20.1925868504071]),
array([19.753843601326082, 19.711945944471964]),
array([19.06311251202338]),
array([16.713906108836913, 18.47560666422467, 17.119786325653294, 16.76393837339465, 16.33725683641064]),
array([18.939943999784322, 18.429392097644858, 17.745826166572545, 17.653084930369033, 17.289767367881126, 18.25433860502469]),
array([16.705508849058546, 17.43560486158198, 15.290852666408785]),
array([16.73608028008953, 17.387903210998243]),
array([16.02682061868148]),
array([16.21861837650919, 16.07443860299396, 15.825398188165474, 13.859489558140893, 14.44121080227956, 13.588695740601498]),
array([17.34755513076856, 17.315093130440808, 16.887414113736785]),
array([17.065740066955197, 16.890139365199026, 17.08189035864221, 16.538138950072113, 13.954937705293325, 15.099383699896139, 13.565523502772315, 13.638751257997942, 11.151053780282187, 10.704368316291768, 10.080209124996584, 11.041447253076674, 7.2932581604964035, 6.052511309301208]),
array([17.131047512044763, 14.297955092922802]),
array([16.09564378909173, 16.329973249472104, 14.70737112709783, 15.568691075968022, 13.15552175038852, 11.860628208839188, 9.25267237091776, 10.039798321367488, 10.074666376569079]),
array([15.974723248343029, 15.017275494641131]),
array([16.009166339208807, 11.308258878610937, 10.70508755150786]),
array([15.617433423055083, 15.820099983648724, 14.785909569687515, 14.895541293763635, 13.937660175586915, 13.802098977619695, 13.178473707144914]),
array([14.438339693961387, 13.037027774043343]),
array([15.491367019462253, 14.446786736388658, 12.865028711816091, 11.256002261744095, 9.45104665542163, 10.747703345963858, 8.554043383248288, 11.447589717127581, 11.002024813974504, 9.323367480030146, 9.477464971276714]),
array([12.92076974445492, 10.280646354387057, 7.750211818989112, 11.241364780193331, 9.75879112860216, 10.6282820822314, 11.00575530878532, 10.826255766689094, 6.845798210709774, 5.231967145238614, 4.890949797832816, 5.152613998930499]),
array([14.545944125689482]),
array([14.652160513208678, 13.983618355399082, 14.13746777249432, 14.251483955002387]),
array([12.796948063412929, 11.784299571993191]),
array([11.10158447687542]),
array([13.250234005848247]),
array([12.686565287929307, 7.761860919252502, 7.809991894219172, 10.81906117164943, 7.592839056257929, 9.59331849853738, 7.4694711790270905, 7.857214433427991, 9.26836530971785, 7.7518983651863795, 8.596747694194011, 9.153673785245505, 8.995011500490879, 10.591297187149504, 7.652757520653223, 10.985233731748131, 6.197745421237097, 6.8261316273063635, 5.381512870415968, 6.987277522539911, 6.747731230051389, 4.777012908072701, 4.9339635781904, 4.80178908337438]),
array([12.650128607648398, 11.308490462965388, 8.883614720916329, 9.409060873364155, 10.332977410118783, 11.574152552543653, 11.456070647205935, 11.207938506540906]),
array([11.78674471765196, 10.710861916637597, 10.215509289788166, 11.524014512038152, 11.494872250209701]),
array([12.212686898522664, 11.904481424364333, 9.926619482984547, 10.522104315713099, 11.298597555652107, 10.526459133075026]),
array([9.820578768277105, 11.359404240556874, 10.960291049734193, 9.939372947074478]),
array([11.77263702870461, 10.061861046644871, 7.887479900446665, 9.91494673164829, 10.38871664586965, 8.09178058156492, 6.184195271873888, 5.735267374929625, 6.137225599605772]),
array([10.45313289678212]),
array([10.121779694807362]),
array([10.57706323321073, 10.951583753142826]),
array([10.949880754446175]),
array([10.418210275958053, 10.645356807880258, 10.080254400847515]),
array([9.292918306534872, 11.329678818377182, 11.355632459665745, 10.831312229761334, 9.321965859788103, 10.7599080578194, 8.795043216636406]),
array([9.536049582553812]),
array([7.9302456293714805, 10.912366658491084]),
array([7.289988756867169, 8.604538313676118, 7.802649085162539, 8.284549288805303, 8.971606245178732, 7.337965044510879, 10.740398006740458]),
array([7.847933996652612, 10.508850864150423, 10.072209698196065, 6.305814650106713, 6.201386925856409, 6.125923046617704, 5.136732239735922]),
array([8.658046359128516, 9.43451717102457, 8.271008220281583]),
array([10.383196030096606]),
array([9.86592954815734, 9.378354502701583, 9.555911139218287]),
array([9.33420785413638, 10.023295359562027, 9.105929329338025, 8.63813176722486, 8.417977683934446, 9.532942926252252, 9.374409886738597, 9.311935900507335]),
array([7.471112235538218, 7.559205401629619, 10.209042780913295, 5.683212520399518, 6.765684333519767, 6.30870261484294, 5.1391854795553655, 4.407321506971247, 2.0544756603293743, 2.497424225371693, 1.4499022301609814, 1.7879998524197769, 0.9517724456552416, 0.5375378704090339, 0.0]),
array([9.498641740176543]),
array([9.169240317335944, 9.122203377958854, 9.608496965178066]),
array([7.517355679172951, 9.24596729446755, 8.902124464132017, 7.9528082475698, 9.214704973074928, 7.447965084277545, 9.283133677842269, 8.817761218431624, 7.741295772697054, 7.144356952130431]),
array([8.93396205440918, 9.227402894031444, 8.605769030568856, 7.36888462179247, 7.822734208494661, 8.095925284151118, 8.105058084738506, 7.88944667013783, 7.709842057959956, 6.734665476179336, 5.951399373821796, 5.494298508005687, 4.602883161935796, 4.321910900028384]),
array([7.980349499723259, 8.141127817437638, 6.523951550092719, 6.52927104711768]),
array([7.587066839176339, 7.796602740719292, 7.924867367780274, 9.047863983926373, 7.458187542280131, 8.23296543760269, 8.961787822221645]),
array([5.834119877035447, 4.977984381063533, 4.548345130888476, 2.6188323285599773, 1.2920774716501966, 0.0]),
array([7.483899368832571, 8.02727868528046, 8.065131554052266, 6.58862671295435]),
array([8.670950741605358, 8.74773881976014]),
array([8.433175717463161, 8.230925073812886, 4.64096196479945, 2.3899607438877437, 2.5523837653265966, 1.2326368714927265]),
array([8.221520665200025, 3.2435021344945745, 3.2461167525595185, 0.0]),
array([7.294040764714427, 5.443760844673222, 5.899947069982728, 6.868269658884271, 6.027738590585077, 3.75272100335269, 2.890919479213191, 2.703031565516277, 1.8142590338913767, 1.6205342563689378, 1.5402975315571228, 0.3675216789897693, 0.0]),
array([7.33885754808726, 7.845378582151348, 8.571722908653632, 7.2557726836666445, 6.928846368969402, 6.46030539866004, 7.139251218789199, 6.466062044212384, 5.117559177430506, 4.81549843987349, 4.553247465170177, 4.232532678588408, 3.064930023384208, 3.2739982545953725, 2.0620950845204695, 0.8170266273537554, 0.8101821813407946, 1.075700942680319, 0.8972185950397895, 0.5475472447537935, 0.6579224643621308, 0.19827724616908138, 0.6947776689303029, 0.144114757660469, 0.6795155805321011, 0.020165431971111356, 0.0]),
array([7.743224386520731]),
array([7.780153054539851, 7.001855177109328, 5.6791897656728345, 5.954253988848085, 6.660695987316758]),
array([7.368170939473977, 7.145274953636111, 4.66741975575226]),
array([8.127148366478389, 6.2491193397062]),
array([5.788027448271371, 5.284133560474635, 1.4299146102109348, 1.1710047683789648, 0.6502113756864029, 0.7581509999585875, 0.43531392400334373, 0.0]),
array([7.490336997565524]),
array([7.034636491626204, 1.3034965273211498, 0.0]),
array([7.375700945896627, 7.956023547011962, 6.782861996404057, 6.116428667822722, 7.05655154492158]),
array([5.967835788007331, 6.371423206169553, 6.603419668791817]),
array([7.582825727383014, 5.652687823691981, 3.9947878856549695, 3.898783591571333, 4.2003987643100436, 4.606854567304624]),
array([6.768257011746419, 6.628579777306105]),
array([7.3009183120773455, 7.843725011374969, 7.795975051576713, 7.324305237548296]),
array([6.739408162189965, 6.000636109112207, 4.459793718195936]),
array([3.58381496120864]),
array([7.238752286184623]),
array([5.071233980551151, 4.966834277271586, 5.012210518241674]),
array([6.335127785477276, 6.032090325370438, 5.8715007638878625, 5.997011762123906, 6.506925152773564, 3.6457156155378767, 4.934197731344152, 4.119173003441748, 3.109838495179766, 2.827725292072648]),
array([6.732172709922514, 5.439254799810707, 5.925306461708235, 5.952142985390959]),
array([6.289420021813266, 4.247434481883911, 4.888435917168968, 3.8385436081082087, 3.518978046055827, 3.325487553380696, 3.3878495677920575, 3.186942484169845, 2.9136727294029643, 2.079030656694457, 2.1922809149483764, 1.4468712153041479, 0.0]),
array([6.733983919996248, 6.627394562675561]),
array([6.74096447112277, 1.4720001591595442, 0.0]),
array([4.571556763102573]),
array([5.820729972096655, 6.488674205055116, 3.6176000185394446, 4.618959032200138, 4.284223715230416, 4.108729088709687, 5.049003164564544, 5.2872528214565975, 3.594031393703473]),
array([6.204014050726896, 6.171088896116772, 5.797830310597734, 6.1981105827103615, 5.388926129163884, 4.750462469790497, 4.529852844187793, 4.81458615178166, 5.08380925302912, 4.169671042449853, 4.480092886195841, 4.384016370314428, 3.3603682813953877, 3.321914011333546]),
array([5.671760210426484, 5.757779604393798, 4.889201074733917, 5.163272608742619, 5.133229234859144, 3.5156350368583564]),
array([6.215701301778673, 4.893531113649371, 3.6880336032473355, 4.784794738611637, 4.376529776297275, 4.8580033956225055, 3.8664792635399587]),
array([4.03946747815117]),
array([5.698251104542415, 4.077230103916673, 5.046450002279013]),
array([5.553848616482245, 4.26237819349144, 4.667808258485093, 3.884853672819913, 4.092337532415471, 4.063425123719792, 4.124724095179365, 2.685736205421495, 1.908169943210394, 1.70139936163942, 1.6624311642199818]),
array([3.958241140168403, 4.938651744239971, 3.1143020710543086, 1.2145848870611764, 0.2540980095025095, 0.0]),
array([5.449648025479272, 4.058808260168835, 5.1971542973377485]),
array([5.400818452869187, 5.234164706326629, 5.022363979009438]),
array([5.694050785333238]),
array([5.3035374199488725, 5.241774704962925]),
array([5.654218815947634, 5.189766545351902, 4.872199698951941, 5.055351761232147, 5.041183767963975, 4.3041597463858405, 4.53238974107372]),
array([4.3652662808164475, 3.1984054733499594, 3.23692691972653, 1.9409965697659435, 2.534337982828098]),
array([5.117225966150755, 4.669890898249881]),
array([1.6738233363324917, 0.0]),
array([4.935326386927935, 4.296822052696543, 3.8389873955967975, 4.094250515772346, 4.139706454491556]),
array([3.991466479737951]),
array([2.70295052670508, 2.9744095839204814, 1.1598630993810692, 1.3154170174806041, 0.0]),
array([4.792251942078855]),
array([4.11170385045373, 4.553817324215418, 4.65597369198288, 4.07280106164016, 4.425782779571515, 2.9431888385340175, 2.4455819470929225, 2.4083035281407654, 2.2643825525165013, 2.095792359035813]),
array([4.658613284165797]),
array([3.9213511656336832]),
array([4.298428110120773, 4.35641073635461, 3.322005723532646, 3.284104874010173, 0.3878554932399353, 0.6399774246064496, 0.0]),
array([4.369947837667184]),
array([3.1783401188880815, 3.001202448052406, 2.6212889656792124, 3.4173772923450856, 2.2998158465443908, 2.239789350122516]),
array([3.674663388865566, 4.153365337672273, 2.922316764189117, 3.475684397815302, 1.9555548122216466, 2.1993623537372877, 1.1519372573543087, 1.6971859093492048, 1.0970743767782487, 0.4053996843815854, 0.0]),
array([4.062037706763549, 4.150311323042921, 2.855988312801871, 3.013289397123419]),
array([2.9326181971283427, 0.7134520642771874, 0.12913674253261165, 0.44359575107588933, 0.0]),
array([2.474639493047136]),
array([3.9302204914035763, 3.282702122053021]),
array([3.9552073937762113, 3.6997083364257035, 2.1178242168929904, 2.420120017681683, 0.4027999076281688, 0.5583633465894327, 0.0]),
array([2.789575754740565]),
array([3.836714748460984, 3.592974555251023, 2.3574129964395856]),
array([3.6333677801149804, 3.1355091080043085, 2.6471856013622648, 1.6132854389853837, 0.2780774911674546, 0.03693827780445404, 0.0]),
array([3.7563095944338345, 3.9140019054944135, 2.8830462098901615, 3.220866962445708, 2.945109456702558, 3.1856483823856436, 2.2026391703815564, 1.2760729837489653]),
array([2.6213138542336765, 0.846055561970572, 0.10048136813393098, 0.0]),
array([2.8542498144651063]),
array([3.406064407106102, 1.8256286957096162, 2.348456832730978, 1.4791546118786862, 0.9779501572501983, 0.415322779758756, 0.2722132865456228, 0.7266678044007353, 0.0]),
array([3.0609122929749515, 1.9746767822002564, 1.2038429320436832, 0.1487150944649619, 0.47803099004414923, 0.0]),
array([2.4129133514555856, 2.4670563626569977, 1.3386346066777213, 1.7741631657579524, 0.4582091776498846, 0.5863124177600569, 0.0]),
array([2.931795925275208, 3.1450874341805863]),
array([2.940828657050406, 1.9464630447590858, 0.061109214403466436, 0.0]),
array([2.5882702862300304, 2.450499241568654, 2.2175118316680793, 2.272673851238579, 2.3974761951605155, 1.6070350901270805, 0.9434914585298577, 1.1436701498587707, 1.4760190659663228, 1.424448129001388, 0.3215638504111263, 0.588597381277235, 0.020654226369240566, 0.0]),
array([2.9052652372639307, 2.8018099437482125, 2.1552845361032644, 1.8351512311295357, 2.405134431987257, 1.4456172826230005, 1.7313974200284947, 1.4103856146648754, 1.4623872545834726, 0.23768857501313667, 0.5523853569037287, 0.26229233778483463, 0.7585343587146007, 0.7551655858722267, 0.24497185717913705, 0.0]),
array([2.897353386688342, 1.8218551728141592, 2.0954346978074256]),
array([2.320420240077213, 1.5183402760738143, 0.887482669205553, 0.07164481625789546, 0.0]),
array([1.8500188376032143, 1.9844022176540896, 2.006500003007651, 2.2636369549765467, 1.6116070174687516, 1.7057811800649338]),
array([2.1652686162935533, 1.8844569612863578, 1.1687952553765943]),
array([2.453732966844219, 0.9534284902242919, 0.8884754547018885, 1.0515365420048446, 0.7784340123194402, 0.6762206814076807]),
array([2.055192177028255, 0.7860146295027182, 1.305517775849566, 0.6698308709492786, 0.658558154874635, 0.05301954549317099, 0.0]),
array([2.3451002233348386, 2.0442905739746746, 0.21171981870676793, 0.3202595367599722, 0.33133508729817557, 0.0936380470034027, 0.0]),
array([2.243392707961456, 1.2312315360523183, 1.1213065743856871, 1.7544417415561464]),
array([1.3643062967294843, 1.5951179003140812, 0.764537656446007, 0.16382170521117212, 0.4824740124278543, 0.4419563733990811, 0.7738567706176215, 0.0]),
array([1.7709046356653908]),
array([1.3303061181165354, 0.24479055262751426, 0.0]),
array([1.830982373275875, 1.6500475837754192, 1.416684959012935, 1.0628533117529173, 1.3278196062990748, 0.32300832794564854, 0.6414070636658951, 0.34495135361799034, 0.1974932292389986, 0.4254786632138644, 0.5705215292757481, 0.05168834635990775, 0.0]),
array([1.872035617164849, 2.1573892750129544, 1.8628796717930025, 0.9812001152646601, 1.7743964829444903, 1.2534179687899747, 0.9337304329179102, 1.5303919069285459, 1.1399119276927965, 0.5183686811942215, 0.3873942502322269, 0.0]),
array([1.2611537499303067, 1.5612991565516974, 1.092297876258058, 0.738157784449686, 0.6683339822323149]),
array([2.013915405847186, 1.4038992262921854]),
array([1.9075702067083604, 1.6546140407276635, 0.6720047490247114, 0.0]),
array([1.9081647773176824, 1.3508132015911452, 1.3746502539040155, 1.441808358290841, 1.5466422341059127, 0.8796629111374791, 1.2266141555886696]),
array([1.3368718171229081, 1.540262133984282, 1.4505097600788888, 1.019276539837795, 0.6548824192658229, 0.09051006495997632, 0.0]),
array([0.8924232548183777, 0.0]),
array([0.8382454629316077, 1.5247385793574122, 0.0]),
array([0.06860697511664929, 0.08241298461549096, 0.0]),
array([1.209031039869573, 0.6162636348148941, 0.6477193204698053, 0.2361234334841973, 0.2173387565996202, 0.5902577067028627, 0.022244649159205584, 0.0]),
array([1.0886267081908731, 1.2104767740677806, 1.333055746534598, 0.3835598478315422, 0.0]),
array([0.8865868480592332, 0.0]),
array([1.159015786376744, 0.7972980700755854, 1.244607018879058, 1.1247718574810004]),
array([0.8277687555463671, 0.728196925381526, 0.4426840731944207, 0.0]),
array([0.9561818847592258, 1.059660761255055, 1.278555948058286, 0.39051044257692186, 0.3974059193782109]),
array([0.7312005298762609, 0.6405480943941706]),
array([1.1663407210445176, 0.8611395175477561, 0.25417501957985633, 0.030954961384158014, 0.10108197820281142, 0.007804940042437869, 0.0]),
array([0.8291414459804463, 1.1048649043300698, 1.0450458963165363, 0.7030764045348961, 0.6133932683626989, 0.20419981340389648, 0.7204963055535802, 0.615114971274548, 0.0]),
array([0.5664749143299896, 0.4077325617048133, 0.3849663939435115, 0.1965779974703057]),
array([0.5752382291142968]),
array([0.25938822967821995, 0.2709571383253472, 0.0]),
array([0.4539831990242576, 0.3724657415541458, 0.5088000527319669, 0.7074990776029326, 0.0]),
array([0.17109798321206582, 0.0]),
array([0.45442117027871365, 0.18132457062236762, 0.3833314826255745, 0.5876106766187594, 0.507429240669735, 0.707900874045161, 0.4726968277220347, 0.03392605858803083, 0.10418640050411787, 0.11226000926928695, 0.0902438118607155, 0.0]),
array([0.5717699298074124, 0.2970820544606608, 0.19945286497922643, 0.4557972578497958, 0.43137198553445005, 0.0]),
array([0.5668790661840645, 0.17745958443161158, 0.5556478405824203, 0.04033049237172476, 0.0]),
array([0.3273848389670272, 0.1735717771936231, 0.4928783277332478, 0.0]),
array([0.138962598320248, 0.3933032076660322, 0.14860265632566627, 0.3389315898914904, 0.08931266959041032, 0.0]),
array([0.13054987183447095, 0.24260258965249487, 0.0]),
array([0.22554226876776381]),
array([0.18377419858140687, 0.07484299389532564, 0.0])
]
d = [data_1]
names = ["8"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T7', 'T8', 'T9', 'T10', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'T35', 'T37', 'T38', 'T39', 'T40', 'T42', 'T43', 'T44', 'T45', 'T46', 'T47', 'T50', 'T51', 'T53', 'T54', 'T55', 'T56', 'T60', 'T61', 'T62', 'T64', 'T66', 'T69', 'T70', 'T71', 'T72', 'T73', 'T78', 'T79', 'T80', 'T81', 'T85', 'T86', 'T88', 'T89', 'T91', 'T92', 'T93', 'T94', 'T95', 'T97', 'T98', 'T99', 'T101', 'T102', 'T104', 'T105', 'T107', 'T108', 'T110', 'T112', 'T113', 'T114', 'T115', 'T116', 'T118', 'T119', 'T122', 'T125', 'T128', 'T129', 'T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T136', 'T137', 'T138', 'T139', 'T140', 'T141', 'T143', 'T144', 'T145', 'T146', 'T148', 'T151', 'T152', 'T155', 'T156', 'T157', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'T165', 'T167', 'T168', 'T169', 'T171', 'T173', 'T174', 'T176', 'T177', 'T179', 'T180', 'T181', 'T183', 'T184', 'T187', 'T189', 'T190', 'T191', 'T193', 'T194', 'T196', 'T197', 'T198', 'T199', 'T200', 'T201', 'T202', 'T203', 'T204', 'T206', 'T208', 'T209', 'T210', 'T211', 'T212', 'T213', 'T214', 'T216', 'T217', 'T219', 'T220', 'T224', 'T225', 'T226', 'T227', 'T228', 'T229', 'T230', 'T231', 'T232', 'T233', 'T234', 'T235', 'T236', 'T237', 'T239', 'T241', 'T243', 'T244', 'T245', 'T246', 'T247', 'T248', 'T249', 'T250', 'T253', 'T254', 'T255', 'T256', 'T257', 'T258', 'T259', 'T261', 'T262', 'T263', 'T265', 'T266', 'T267']
def get_taxa_names(): return taxa_names