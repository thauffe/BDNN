#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.70196399790605, 25.882002145339694, 23.5159154189167, 26.33434478219187, 25.49407328995976]),
array([30.769253823256232, 30.964566556393887]),
array([25.404233636170467]),
array([26.089528093631163, 24.946020429844683, 23.727815996654, 23.285930438055757, 27.35433745168978, 23.402200740431226, 25.60497039408674, 26.790700951324283, 25.772181715872016, 23.131234843302632, 25.602426263456152, 26.433229335842334]),
array([26.738534325539216, 27.50931458719474, 27.83130787690702, 27.727395745055126, 27.729030155642768, 27.682270720715902, 27.707898741607156, 26.7238855001857, 27.623515034577675, 27.302492124685227, 26.222294716584216, 27.05103518508477, 26.671586832498786]),
array([27.807313699174745, 27.791937455517324]),
array([26.70739415120056, 26.323194810783505, 25.994228982904122, 26.298968547169757, 26.563892398728832, 26.468092386698196, 25.44766369590182, 26.210998844782623, 25.73719041685577, 25.31974818743543, 24.7588155141119, 26.712063580476702, 24.43176880327287, 24.946298831527276, 26.58472290415338, 24.209264127679667, 24.06738279840456, 26.64676612043541, 25.999367055558388, 24.411396160620814, 24.585220928849022, 26.276436941930868, 26.20195202998987, 25.718583691625852, 25.376431044426933, 24.97523568540911, 25.72325908129068, 25.442395851600555, 24.58758718559873, 26.719950391328204, 26.072689166967937, 24.68995124234634, 24.1443794869384, 25.294967673652998, 25.233659544204052, 26.50848396598067, 25.487762899529635, 25.232148200634928, 26.102588489095808, 26.022270307151857, 25.69814751068613, 26.769786406666817, 25.863762503250463, 26.195655979010365, 25.847507943466255, 26.470713604542244, 24.158418157780314, 25.049222639337852, 24.79269304256653, 24.462486284177984, 25.334673403212825, 26.569239135737728, 25.218680374544356, 25.06947664593736, 24.593681577447008, 25.53502136131096, 26.447259091613308, 25.978672282603824, 25.740753110804093, 24.892248276005194, 25.117665085993263]),
array([25.3029785575806, 23.077601171760783, 24.421614887670188, 25.72258942517194, 26.200707871649975, 26.020407238679308, 25.691769218802687, 25.992794752393767, 25.90141840509889, 25.553614054510717, 25.947991208609515, 23.39441805106872, 26.2767748778144, 24.708474794817, 25.111299200726496, 23.09220348231262, 26.164502675346775, 25.15556973536853, 23.13776937053795, 23.50475984559425, 23.675735589549042, 25.649603020531362, 25.956766816513067, 24.19968943366676, 23.27930275837962, 24.594626987660828, 24.542207206171256, 21.487952805141088, 22.856271923563973, 21.08760238166613, 22.486344513882432, 20.70403939034855, 22.4986883311516, 22.536118803481106, 22.935336360273435, 21.8073466121651, 21.567061464849825, 22.859814623880567, 21.021383472759922, 21.85598191800261, 21.671829202441568, 21.88842448542999, 21.992023672332905, 18.99346232040155, 18.255764938891144, 16.996278842399164, 17.385669666568923, 18.446814902495937, 19.216888801720955, 15.982853484295031, 18.534969579005203, 18.422300001246352, 19.70860151927324, 16.887201317087815, 17.339240450756748, 15.751945218834186, 14.9253757177494, 15.417930379400909, 15.357523260781866, 15.309326575408484, 14.243494319659774, 12.84062643629441, 12.504698409940248, 8.052888619561582, 7.853753935942697, 7.9188811391818215, 6.623100706432865, 5.930997091149438, 7.112567560120831, 6.910927166370862, 5.606407411936513, 6.810723968742946, 5.6422292040345265, 6.71090701452135, 6.014238260327517, 5.916007501155164, 6.263753492974079, 6.6610247434490635, 6.83857337487298, 4.74639220408947, 5.29656319851255, 5.114842669729118, 3.036238782767084, 3.1955152443702377, 3.4656506708493247, 2.7173983259129093, 1.9337266336332566, 1.9643515783106453, 2.3399192512890186, 0.7834167649472024, 1.5539150815232516, 1.758882507429966, 1.3809958431635168, 1.6674477974443032, 1.777148018943866, 1.4327250631428257, 1.001370722154334, 1.223061875641653, 1.4812792289946493, 1.1332254964029187, 1.563573474265322, 1.7853453064525453, 0.5660424387675642, 0.6440661102220159, 0.0]),
array([25.138826969178233, 24.507354116182402, 25.258713564360797, 25.082453785038528, 24.951186994184823, 24.845182611497055, 25.187291445184545, 25.351860496086427]),
array([23.851591710448155, 24.53098953266306, 23.913587352798288, 24.814294670771435, 24.124794385086336, 24.016236297316155, 25.03065233730468, 24.42976497301083, 24.856385031658764, 25.00627881230246, 23.051183909657215, 24.459833387103334, 23.431164919246722, 24.60270821844846, 23.70202925610751, 23.535919215632763, 23.88320225731745, 24.144974242159165, 24.271015671149264, 24.403887957594282, 23.17869350563709, 23.387393049567894, 23.10379321849366, 24.8777328584455, 23.59786962732031, 23.160101636107612, 24.757598230390457, 21.580531541911718, 22.627927763039455, 22.55880614201725, 20.797975241574992, 20.579560953073894, 22.063128870835442, 22.328771296558198, 21.362041183010952, 21.97137849685593, 22.964249338839707, 22.62556477282157, 21.17563192728524, 21.323386990452708, 22.959821379415978, 22.455092100593664, 22.431516276836895, 21.289797794079295, 22.9789497072162, 22.760358772659956, 20.95821733888344, 22.567569325279276, 21.339387580663647, 19.45315983083424, 20.071756726097185, 17.398168119726034, 17.514334872163325, 19.218529914318005, 16.189639876740667, 19.764771223739377, 18.50746343933957, 16.222036836850833, 18.563818893340233, 17.58578901394154, 18.66275304110631, 19.214258045501445, 16.04161991361957, 18.61027699177176, 16.44337655955427, 19.412102911666846, 19.48742197300788, 20.300830649506235, 17.527049290569796, 18.792687243554788, 16.576318623216824, 18.32475171708874, 16.112175377221195, 19.61070094480705, 16.55725443923842, 16.968612068020917, 16.993979354441088, 16.663818970317514, 17.038222760484796, 17.640795162307576, 19.85333539675973, 20.294607977674104, 16.32381012026388, 17.683573094092303, 14.810803255148203, 12.69066205334467, 11.806531578660426, 7.601426076524194, 10.524175367559984, 8.319896184376317, 9.61099484395017, 7.93867834412997, 8.436296871323488, 7.182409196658392, 6.663638695222773, 5.994540093342062, 5.448207967852707, 5.576772211983992, 5.77821988815836, 5.73533489032668, 5.971888924939993, 6.7472903809035385, 5.50781432116763, 5.960876890388422, 6.172195751041477, 4.614398172390761, 4.4773442689059175, 4.630578737249077, 4.411477816347078, 4.300422202413171, 5.312831789448844, 4.527115374145438, 3.6631463132644724, 4.432074155375881, 4.157851877296011, 4.790475525803792, 3.9909949157957207, 3.992845836354321, 4.42815783755834, 3.4395813276566902, 3.284016383606719, 3.3041535573983145, 2.711157685125329, 2.913842945906924, 2.7873849748648003, 3.344522743231651, 2.962628418400861, 3.0168190971933764, 2.7505130156759634, 2.6757911983091143, 3.5876000930930654, 1.0211238992212803, 1.1644057965235746, 1.0232136256335216, 1.7603221609952688, 0.9469750992405653, 1.1687028021630068, 1.2259001376238197, 1.1635182849712922, 0.919142521441137, 1.2019116998164736, 1.2080115628939452, 0.7904979389499847, 1.028891959757129, 1.3493026122452931, 1.0791624569664466, 1.4033741493460674, 0.725529136071674, 0.47531600582071254, 0.4652034752075036, 0.23020479468609578, 0.0]),
array([24.54227604686582, 23.449384548150103, 23.719353556648542, 23.2862373816516, 24.42356874986541, 23.36166695478698, 24.323363741522652, 24.384877372557558, 24.38388940735179, 23.33038895322516, 22.41310659708898, 22.772446792753296, 22.400184025311415, 21.159328646493844, 20.81374614760472, 20.879730550771495, 22.001294658288202, 21.29965949628203, 19.970042998324057]),
array([24.250793500824052, 24.327544670624675, 24.416231538573307, 24.13631170796092]),
array([23.91095459423718, 23.295930812766215, 23.843254382573186, 23.933430577333883, 23.809114310779073, 24.043002082461996, 23.395966398744438, 23.751055434365906, 24.143291592037116, 23.175616868751877, 23.923357729008845, 23.691848351384717, 22.499090674873024, 22.270747378623977, 20.7084902300303, 22.996325018107196, 22.808532905086928, 22.49157459917179, 22.96128159196681, 22.560082223117426, 19.727346919761427, 20.266383880335656, 20.254725027540225, 19.425363374901124, 19.803147275548795, 20.143223404483503, 19.712440070822673, 20.344970773128804]),
array([23.148780419075027, 23.50721774518923, 24.15153948830038, 23.106413203010888, 23.5028107725993, 23.874669869385503, 24.1231244677175, 23.291588009626725, 23.51250759997079, 23.62250441319007, 23.164756419525286, 23.29848154174485, 21.519922357088497, 21.93851041396543, 20.717980268310637, 20.912821175097438, 22.1446930591025, 22.358093307412954, 22.04659666284678, 21.50113609888665, 22.335715170268752, 21.413454328670745, 21.969623820824687, 21.60019780159337, 21.90043228953953, 22.267333682786436, 20.573505934968015, 21.05226099312622, 22.66710687476489, 20.780667942424408, 21.76755319653376, 20.713109732344304, 20.958374887096454, 21.928697023424686, 20.311966291578834]),
array([22.984387304110165]),
array([23.34388029449888]),
array([21.975777545291614, 21.40256814921218, 21.20969573836509, 21.998435860590273, 22.247352244447068, 20.786646502721474, 20.472755420522542, 21.43605061716436, 22.17947349971856, 21.968154614309324, 22.50137957925032, 20.619972506222954, 22.2927227259307, 21.329315401420068, 21.96448403655687]),
array([21.353935094842623, 21.154588203878678, 20.789029302400237, 22.066836392890764, 22.467458551250115, 22.209384458358944, 20.644961784211386, 20.698525340996845, 21.26530314055087, 21.25894043611101, 21.126674817180177, 22.13048887882193, 20.615455510853835, 21.63152042397125, 21.70456601708853, 21.426218704170182, 22.393719768111136, 22.149119094050036, 20.97703108632594, 22.32026273075349, 21.253538175593903, 21.987034325153292, 21.691943991213947, 21.898177839620065, 21.760870931246853, 21.015955894372297, 21.88090048297249, 21.32887048739136, 21.514096675491846, 21.129474111006573, 20.569499034384247, 20.98822720071899, 22.04937222388895, 22.426074573124364]),
array([22.344113258736083, 22.184362593263877]),
array([22.356119856761815, 21.457340180653503, 21.254377445187885, 21.78482734054974, 20.90275734843344, 22.24467152866978, 22.308221680893666, 20.989058712468722, 21.45910754126348, 20.05533943165781, 19.29949048608882, 19.99618239873929]),
array([21.035151422823027, 21.519479783266846, 21.879372116513306, 20.78551457324228, 20.82311567116925, 21.023129740668047, 21.830202787890773, 21.222659426280046, 19.89373582110184, 20.426055352110282, 20.35306373464653]),
array([21.48161254798776]),
array([21.618491954739294, 20.883590917838603]),
array([21.540967791100503]),
array([20.49772466956787, 21.082633303368773, 18.68555380286212, 18.294673408814457, 17.75790814990671, 19.136268273169417]),
array([21.451150303371037, 21.183574289466232]),
array([21.079879790267977, 21.157276312474742]),
array([20.63343885213855, 18.403703073649233, 17.91984137623563, 18.155036602542488, 20.323916841729314, 19.032114644534563, 20.13255457723475, 19.730808103522797, 17.549078601727263]),
array([20.68226380963539, 20.658483639192326, 20.547071087388623, 20.972735713709977, 20.17923315649023, 19.52703467918402, 19.30877500619784, 19.15448910390429, 19.41584290921895, 20.422119349222896, 19.842575805913494, 19.49258769260598]),
array([8.226054394545208]),
array([20.65359551236133]),
array([20.681425073035925, 20.44382291752104, 20.05075140878839, 18.88672269474333, 20.292910228054247, 20.09189725760326, 18.68162988642693, 19.852069209589246, 17.9257940422013, 17.385670186721914, 16.893400922946764, 19.553334417995494, 19.698052160592322, 19.97676446731861, 17.62061399452169, 17.31031253166303, 20.300073378422407, 16.292579749140454, 20.20976721342398, 18.04205723461626]),
array([20.448391099047008, 17.270499478054404, 16.546768016854998, 17.248031413368047, 16.85627546867803, 18.922657639333657]),
array([17.336257572313116, 20.09876817759362, 17.035489373212418, 18.34105266566195, 17.326960384941025, 18.69928325021119, 17.980147866061067, 18.420725994880215, 17.64922339819995, 20.093071267564984, 20.17307616419548, 19.470891840121848, 19.800112002774412]),
array([16.532870142660236, 19.397650595498185, 16.075951526954015, 18.881027474124842, 19.560502350126608, 15.300688446814593, 14.657782365881525]),
array([20.211545223928628]),
array([18.593946861452732, 19.009673094386372]),
array([19.113667421368522, 18.82858890580626, 19.68248875008354]),
array([17.767845419638576, 17.813942633005887, 18.508303031217554, 19.47297915501163, 18.444311561197036, 19.745495404784563, 19.834403645279572, 19.348558618586, 18.909986872577864, 18.041896838426066, 18.3438581907001, 18.31060804463204, 18.98944426428112, 19.1596991159495, 19.10597704721295, 17.786429933604836, 19.814166330377287, 18.04455259829796, 18.66163462565508, 18.961363002934277]),
array([17.703015667816864]),
array([19.400732355187376]),
array([8.753386503606183]),
array([16.42441057528684, 16.32648030819958, 18.046021872524296, 18.566188207926164, 17.43346232034983]),
array([14.430276050420135]),
array([16.26445915872792, 18.62217504494096, 17.236832279226697, 17.207957660489132, 16.866740244684657]),
array([18.70556179184913]),
array([18.178060544500568, 17.934619023961293, 18.853599282823495]),
array([18.591923748133976, 17.300344634486713, 18.6914805632313, 17.50628223566954]),
array([18.23895768652028, 18.523211953692314]),
array([16.312512761627847, 17.439941074242938, 17.25573401208632, 17.668031571090623, 17.048680243697444, 18.18666173764969, 17.583633927872835, 16.05105734581149, 17.93594640334201, 16.879647838965575, 17.75824632980456, 16.570296216930508, 18.273184799096374, 15.618315039053742, 15.436521659032296, 15.812766939595338, 15.736317339100262, 15.091357968385172]),
array([16.049974230541856, 17.417457715333356, 17.99446627213187, 18.42523361957547, 18.170141551518018, 17.242599036326553, 17.72694158245709, 16.423459080411085, 17.478148062296214, 17.96090717764278, 16.98456578549932, 12.731649972436232, 11.647173947536945, 10.852152274644805, 7.0333441769869705, 6.389441578066889, 6.36395355237711, 5.589496124154564, 6.198622010120197, 6.517152548172695, 3.1079654590870813]),
array([16.650276856825766, 17.902462773669797]),
array([17.634761567423364]),
array([17.761791145903608, 16.460292251977908, 16.595323578623443, 16.947186571105345, 16.259842964917294, 16.3314535155129, 16.74124088222204, 16.21501541824534]),
array([17.256537239740435, 17.408690216575682, 17.153005629342477, 15.579907552454403, 14.84089872511582, 11.95537934253824, 9.240601987555824, 7.313207861924557, 5.374450570573577, 7.03385934318678, 6.436397211354315, 7.165480722371588, 5.828343704610561, 5.486965325986428, 6.7094484765812625, 5.601417347252885, 6.686109870435004, 5.500964370290259, 5.801408992708562, 3.6414960857269607, 5.11279812846061, 5.300680539741025, 3.997626293660224, 3.933620920358753, 5.122041493621404, 3.5486825588714686, 3.0566299215734176]),
array([17.459802767378882, 16.71075214575837, 16.07016300617396, 17.487837205743755, 14.44555053729677, 15.740214615346302]),
array([17.314356918057893]),
array([16.698575292383033, 16.057946226687886, 14.841048682938482, 15.459440060434286, 14.053775518277217, 15.00774416302735]),
array([16.56494418867023, 16.3987057164411, 16.80383161262616, 16.71094141516298, 14.599812665382231, 15.034614878334356, 15.896198694146808, 14.893642661490162, 14.352614318248683, 14.947259905751237, 14.577378938974457]),
array([15.971336987704117, 15.973488130356396, 16.06950600190656, 15.695895996361994]),
array([13.363361973557435, 11.461270148030518]),
array([15.56822997172324]),
array([15.31528287535785]),
array([14.535794980514176, 14.903707615555664, 12.666061933841544, 13.67570566662201]),
array([14.109851041811236, 14.0009581466496]),
array([12.029785138404291, 11.440726342265902, 5.703798615571804, 7.21149943946953, 4.8705416989816115, 5.299634411375732, 0.0]),
array([14.718744282769222, 10.907451847842148, 5.401255355846167, 5.853698379193261, 4.822576103349623, 3.306614616375448]),
array([14.85022474538445, 12.986434987286291, 11.052501057175201, 6.662096901240169, 5.656879822130449, 6.706605383836995, 5.573803643391299, 3.992464558120698, 4.433928287057238]),
array([8.447304074624466, 8.044829113907154, 7.167644048273292, 6.756234559687286]),
array([14.078979892012244, 11.897236238225299, 10.693966170094015, 9.412341794099182, 7.9443579251989, 8.417398027744758, 9.252832091123738, 6.210021718748352, 7.061039955105536, 6.782531425056065, 7.222233866492665, 6.347840578170008, 6.66877576038384, 6.099547060085494, 6.67649236043029, 6.671354982352618, 5.822041831230248, 4.727794623818862, 4.686891019812049, 5.327658227512241, 4.811345087312142, 3.287087868487416, 1.5667759037293145, 1.1918079374513733, 0.8541824420218834, 0.8532068291042434, 1.4515985216105054, 0.0]),
array([14.465422628362159, 13.977873256522827, 14.332872954252602, 14.115880589630592, 11.918955848664684, 12.05502337374241, 12.089833344487536, 13.15275268537529, 13.642786770117757, 12.937473408504298, 13.439981312179546, 12.29805711530361, 13.06204159888225, 11.683413660116788, 12.416271150613724, 10.230875973612708, 10.425185503808416]),
array([14.597839630102186, 12.01876756564075, 13.529770966695352, 11.932188283822958, 11.799161752775667]),
array([14.068252845713463, 14.18183333203343, 14.092838414304188, 14.36911062288839, 12.753334397786531, 12.49230606632536, 12.14506653602244, 13.437103406299697, 12.719380292545061, 12.946927969178635, 13.499440957164675, 10.325547589418292, 9.152943583306074, 8.609926835812402, 7.80135233686622, 8.840498704154978, 8.293243181602795, 10.928245952359152, 7.153758326265475, 5.516053189212153, 5.613062416666496, 5.966486226938205, 6.109387471674252, 5.5363484135644265, 6.890907376022984, 7.094695085242373, 6.787098514453326, 5.723588018010435, 7.179571768952059, 6.129342103980035, 5.338741718977451, 5.606486284207275, 6.246910059792483, 7.124064017806794, 5.65308834643698, 6.098357178332465, 7.02578837763387, 7.06616944227845, 5.503155013959489, 5.963812296210642, 5.657086492509096, 5.561876427059175, 6.150484173475415, 6.254572790037341, 6.421301081824983, 6.638730337792216, 6.413314138994709, 5.76324657962102, 6.092714706175008, 5.5359253764465155, 6.396483892449477, 6.082816367174868, 6.03641472827575, 6.112308237021081, 5.899753690768151, 6.791077474593871, 5.200708251672414, 5.135161111236728]),
array([13.833156312752791, 13.58477796973342]),
array([13.111606028077066, 13.178959067950316, 13.005085600243833, 6.264311002222265, 7.046305277819206, 6.50927629623254, 6.311629896521658, 5.996338318710875, 4.73705742111127, 5.31348352474813]),
array([11.683483088195699, 12.283182172952227, 10.116091483903089, 11.599392286967134, 10.820893324369077, 9.247205574429344]),
array([1.2622966619200928]),
array([11.664177838846781, 12.197502121669869, 12.31818386691474, 7.281905173284391, 11.195203211409646, 9.869848919264113, 5.424782503244835, 6.095152556493351, 5.821274056605711, 6.275994230007192, 6.031056675882949, 7.11009387388769, 6.758585731107021, 6.552689473418955, 6.340702731807521, 3.775490153912528, 3.975066744993377, 3.979949306222204, 3.2399823264785965, 2.9662817019693977, 3.2592690350178555, 2.852646968397914, 3.187382066714653, 3.5302797356838425, 1.8672450137788577, 1.5862477424444403, 0.7944125794061279, 0.809070913053385, 1.1439959809954783, 0.987541799725752]),
array([7.198251412204507]),
array([8.941761545170824, 10.9444402676469, 9.764914570344319, 9.166880445340293]),
array([8.929519786059377, 10.805699563446016, 10.071994050881719, 8.19385681681602, 9.20027096204524, 11.607073374227872, 8.293816508036837, 11.475195019941992, 8.36070309650628, 6.977713829087843, 6.715623786271858, 5.384543364233587, 6.611139134390555, 6.32240933265725, 6.87027678239324, 6.075745383818829, 6.708210959071393, 6.849402219711974, 6.256246843897505, 7.235450888808634, 7.231315713444085, 6.826264161760841, 7.175641869448427, 7.039022655085596, 4.537667336923821, 4.64145084318871, 4.805955789118577, 5.216347894008948, 4.935152821155901, 4.903375711488017, 4.492546669401351, 4.03692117658772, 2.8005649606662106, 3.30852650891544, 3.1585777828087522, 2.9008472391136606, 3.260894719708314, 3.5422826307303503, 3.2186005320469304, 2.636662445651811, 2.9454309147554705, 2.9150912683594488, 2.9550722032896073, 1.0713129356519095, 1.1955717747987817, 1.7456383910683133, 1.3365951349388383, 1.1123927580914086, 1.1095952854525253, 1.4197880518275912, 0.6360886538278281, 0.531450715747198, 0.25149424102463247, 0.24152653495942233, 0.5290797452339069, 0.48027376474722405, 0.5823477414937197, 0.08386545144741767, 0.0]),
array([11.728771627431952, 11.073178459545598, 8.995297929671706]),
array([10.016947746291663, 10.196731793973862, 10.977274787342564, 10.2902311192909, 7.081313724421367]),
array([7.786763060145169, 7.262070263842149, 10.015331018257978, 9.663290993895286, 8.464441006919202, 5.7082458379250784, 6.14995479930442, 6.263728637067507, 6.683675414733239, 5.460759873566724, 6.485319180989541, 7.079130258819302, 6.063700555118293, 7.027885638098832, 6.584552458391959, 6.273798725416679, 5.5709896427032035, 6.300364752068333, 5.900472204283331, 5.3437629364777814, 6.994066846275468, 6.719614397104601, 6.48999223915491, 5.031499808109293, 5.293557996848289, 4.726592425377432]),
array([6.4377772099196, 6.458483619791215, 5.5070173191016405, 5.877913453201503, 3.699086927967186, 3.3970280182613743]),
array([7.695981884041096, 6.081253609802212, 6.65324187257491, 5.667617623831999, 3.9716220903532538, 3.6150027302938987, 0.8404158655954668, 1.6294631841507201, 0.17017655463535064, 0.0]),
array([5.855431895277857]),
array([3.1732994586156744, 0.0]),
array([5.672584090633766, 7.2430012584531385, 6.026654974728592, 7.148024609156935, 6.752424565840589, 7.13195067413367, 4.388045564255545, 3.6644107540619535, 5.199203279733021, 2.8340414056694154, 3.496599644495778, 3.3165901683884935, 3.3181545093648293, 1.659899355059643]),
array([8.171417184059253, 8.629670562343616, 8.58169793391023, 7.938164334866261, 9.51268166856387, 6.089127912148973, 5.930493725748166, 7.0266416854657265, 6.510564700085819, 6.371850848841182, 6.538894410306084, 7.01411609907037, 6.443296034000797, 5.806742860399071, 7.109949314754283, 7.092375093019073, 6.566258204842255, 5.44165353653814, 5.965014808878896, 5.945668363015366, 6.552561297554651, 6.10854592404444, 7.226588651402844, 6.0272399750793175, 5.49614220844393, 6.0723092346513265, 5.571993494918483, 5.8588676257269645, 7.0150100741705295, 5.516899463240327, 7.118693594214336, 6.222918322598485, 6.973665270811735, 7.223570133839588, 6.437382357595444, 6.680702004134637, 6.314823053883793, 4.697372305507828, 4.940090406348295, 4.9912087644477126, 5.086410898415671, 4.253224400090826, 3.9783131895571624, 4.329543769002776, 5.075754512412398, 3.9267558780150953, 4.22979835660675, 4.658473328400298, 4.163754588722611, 4.880364385909899, 4.732817479305016, 3.742749491381657, 4.827427205188442, 4.51692029608677, 3.790055832304895, 4.105642064037178, 3.159340693853604, 2.600023884238971, 3.531147402382502, 3.0709249857317933, 3.4170480243475456, 2.609909923682061, 3.4844304251779867, 3.188326977591177, 2.8473573154687015, 2.881354483412006, 3.168816790368189, 3.2465453486882274, 3.1944353738760722, 3.330438211578188, 3.390560402519064, 3.107299508489932, 3.3260955688110405, 1.8368042867380452, 2.3251408172165684, 1.0959892278947114, 1.0896186137962465, 1.1390369030525769, 1.5194813104925546, 1.6408255065083321, 1.1519768069635816, 1.1893111436539128, 1.6172362255152635, 1.081463210663596, 1.3361797404234041, 1.119551957156387, 0.9467229457142164, 1.6270467312801153, 1.7804913795112234, 1.0060944644626932, 1.1787183237694845, 1.087067278555058, 0.9459200804865732, 0.9819356574138339, 0.8208223702795276, 1.2538821546496917, 1.3022229503195564, 1.390476943311346, 1.6433888273417556, 1.216513229693219, 0.7032883334392649, 0.1693941279340372, 0.3330674740478502, 0.4328202399901911, 0.7596083151861552, 0.46948449435658457, 0.5139249912232104, 0.6754016781562947, 0.3510517097070697, 0.1448961556737728, 0.0]),
array([8.034971714108291]),
array([8.531146937179933, 7.6563705214334945, 9.054576664192675, 8.067124679704909, 6.140085017016601, 7.162195937456057, 7.010325686542903, 6.175641168400985, 6.588491351321187, 5.848903597855258, 5.996209201297844, 5.625822040511339, 5.9191837441017645, 6.4222520557204765, 5.956595951492268, 6.718369670126618, 6.659103336340944, 5.73430443359176, 6.390660206574569, 5.716303267566837, 5.670238732051579, 5.60291056104098, 5.489726888128247, 6.420232317496525, 6.635051176096699, 6.130507201525678, 6.534833636693506, 6.094341896410754, 5.894097487174433, 5.927015984569066, 6.528356668469447, 6.92629585354519, 6.421105787506411, 6.972789619939063, 5.540574080587824, 6.583298785226265, 5.222388418735064, 5.075165733397381, 5.316434047964584, 4.2385053029035475, 4.555524974039976, 5.0977045079390635, 3.6727739705717872, 4.413688677574077, 4.911524601184702, 4.463955821902356, 4.154256824368264, 5.1981061529057335, 3.306789101804805, 2.930108815255952, 2.61209192166611, 2.9842900554294376, 3.5214809435937346, 3.5229047138688516, 2.629548314737523, 3.2967633938865197, 3.506091674195875, 2.988974521106611, 2.9638782461439286, 2.9385250834789383, 2.303548539941417, 2.185605360330872, 1.5177242377118276, 0.8611876450549779, 1.2057285704416376, 1.1958958560622357, 1.4383317705663268, 1.747022850334716, 1.7676531030849247, 1.3686214881078917, 1.7138278797189146, 1.1816229523045982, 0.902288507253189, 0.900309801967422, 1.1284327270752552, 1.5847547637290833, 0.9206062288453093, 0.8970707894158408, 1.5996099969329352, 1.0485628185996798, 1.5080365796432618, 0.9821095217123647, 0.7838145799982288, 1.7694035005693487, 0.8855841827590403, 1.1655394862629676, 0.41276632875273606, 0.3128366039292155, 0.736164953569259, 0.5074281177398685, 0.7047851922172552, 0.3796187277235857, 0.07628388826981561, 0.0]),
array([8.270559515364793, 7.952049411863081, 9.091979004490621, 6.950055124320112, 7.1972420920205815, 6.008604307875096, 6.025695159756712, 6.593402416830245, 7.108600020967901, 6.461151037380736, 6.649707279458538, 4.8547623999671465, 3.9339497007706914, 4.534779774203396, 4.38442154428413, 4.581228316177619, 5.254053099227601, 3.327249931585308, 1.8459257662237347, 0.8081231110710472, 1.3261553016405703, 0.9365116987477701, 1.7364204440983555, 1.4583506479702166, 1.7317101480315171, 0.8061253469848271, 0.1662941863697287, 0.3475355852310768, 0.7306940033540539, 0.5486936357347707, 0.4203154864662423, 0.0]),
array([5.948636854267722, 6.085326958048645, 6.547438416773514, 0.0]),
array([8.04151070661112, 6.82851786604744, 5.352923088845886, 3.5864003842069727]),
array([8.041466957441424, 8.042962588729067, 6.376132893788244, 5.927266564245066, 6.9767005141198615, 6.315610135452594, 7.193239644730104, 6.277547994083578, 5.0432528809111306, 4.5592129649765285, 2.615684323156439, 3.3440890933950755, 0.9675035884899983, 0.931254177750356, 0.8442134291688618, 1.7132795476439389, 1.4070760532289133, 0.5168989197327198, 0.23496317313876214, 0.2473904919149571, 0.38025068183238997, 0.5363850459380781, 0.0]),
array([7.7352962194199755, 7.6452385654626696, 5.398942647962483, 6.509859937919207, 7.17128363990647, 7.186688704598538, 6.502617367941799, 6.9393922217365205, 7.2012573834351965, 5.797094715671461, 4.968109144844285, 5.2435362804187555, 2.910453186186402, 2.8792375327237165, 1.727615448064162, 1.5852018295205208, 0.6137546751850751, 0.0]),
array([7.463242254355623, 6.367469914982298, 6.691807016900395, 5.521440621340423, 7.023477752255793, 5.414459523721608, 6.644439769504992, 5.543818563856567, 6.090761448665352, 5.66202869873531, 5.36063077928875, 6.693415181139284, 5.404464409959052, 6.499752689366963, 3.8704860811903448, 3.6325369245460797, 3.957862827985641, 2.977792434578925, 3.280834457902079, 1.566276273198612, 1.3219797519425174, 1.5263351513647383, 1.6971033091527015, 1.6293265802131043, 1.6499262007404194, 0.8084347429273679, 0.7116523215704004, 0.5485420173163071, 0.5103720012289437, 0.3383359840568318, 0.3499053150850894, 0.05922748604679065, 0.0]),
array([6.507359190524057, 6.593818308248368, 6.208663024322027, 6.950658097270019, 6.185088794054051, 6.208093579422433, 5.786196127049069, 5.6487266729275865, 6.402114372801396, 6.666729494242539, 6.749114336740448, 6.873100728873466, 4.749426253240827]),
array([7.135547232506669, 3.323481821160653, 0.0]),
array([6.589869672552689, 6.728205421877419, 5.394344288679964, 6.177429888904916, 6.580617193217548, 4.9113795645515115, 3.6681153323621647]),
array([7.665122906646405, 6.302159680784251, 7.019666163611899, 6.523984863797218, 6.22774299336935, 6.451052125505474]),
array([6.9944864707903935]),
array([7.286938438017866, 5.341786486868179, 5.559316915722656, 6.051554184166513, 6.155072975974747, 6.056814520656097, 6.860286001412218, 4.989318531295569, 3.6840264280740302, 3.608962621438905, 2.990378351643243, 2.7165824537160748, 2.6343315964981366, 1.0262809354500755, 0.8138510668903163, 1.3519044442484411, 1.516342517185753, 0.9865392320438388, 1.4675472575160715, 1.49492834219425, 0.7350791685349767, 0.10566771284306373, 0.0]),
array([6.488916083905566, 7.170395534981037, 6.33485379217524, 6.205564897137117, 6.809631505626251, 5.398920728830876, 5.503484965444244, 5.862373560866485, 6.224031219285566, 5.579139090115036, 6.995670109234245, 5.7313832341693205, 7.006843320439519, 5.7242923705553395, 6.717629339230401, 5.56549444676915, 5.437215348230538, 6.960906922619799, 6.3595135406947065, 5.366428183440295, 4.717918387058398, 4.486216390659033, 4.840195789572291, 5.106794066468562, 4.954183538327566, 4.819537595574244]),
array([6.182377730644308, 5.5551465878510005, 6.008592988483322, 6.3607056519326575, 6.508457133382315, 6.644681899707956, 5.7230329842267995, 6.342208163742004, 5.2525335442513645, 4.670316139505719, 4.2899514930319675]),
array([6.028876089052925, 6.05597456359831, 5.934490988196623, 5.5222246817568905, 5.001113765874151, 3.4919760795081674, 3.4420748151503258, 2.6363690160412587, 1.5630763939402201, 1.6678571211895041, 1.5386641408239154, 1.4379291763051083, 0.5221051483925917, 0.661934841972628, 0.0]),
array([5.355424066469873, 6.111988844682562, 5.55767502281968, 5.640823311055772, 5.892941784078643, 5.388017770432248, 6.385373639502508, 4.513627030222965]),
array([6.596294315876151, 6.298415213217831]),
array([6.036671231041477, 5.522819243699104, 3.7593350867612587, 3.697396244421018, 3.6689654043785764, 4.427874825011854, 3.639257319174699, 2.8487340670985115, 1.896470842148997, 1.1131999018995007, 1.736246425778379, 0.0]),
array([5.450808068129181, 6.38942291921045, 5.461684817825751, 5.953399467432553, 5.883211287425217, 5.751950462464673, 5.413187654395002, 6.438398813294366, 4.663567017721694, 5.187729136111473, 4.910327906308186]),
array([2.6993462986736314, 0.0]),
array([5.7642199212962115, 5.589190704782634, 4.338702977986634, 4.988670364169665, 4.445951046702536]),
array([5.549061503418399]),
array([5.4130162434849804, 4.176781171249349, 5.048529443211181, 4.352402409456731, 3.500710009273682, 3.3120829231563933]),
array([5.5637025150598065, 5.519163191139646, 5.503146808137406, 5.382801953926534, 5.106228195816204]),
array([1.5854934284810156, 0.0]),
array([4.465828251598044, 4.830922662417919, 4.833221148190531, 4.572227992022491]),
array([2.884837613694177, 0.9756049105427816, 1.6502865598681375, 0.0]),
array([3.6797104060936756, 4.4330290780281025, 3.994887216018578, 4.440158059651214, 4.332596411167687, 2.906285876844088, 3.183685275315968, 3.5842972104748654, 3.005415422059652, 1.4508445956935505, 0.9189125586265552, 1.642131107944832, 1.7156911645050645, 1.306255783013469, 1.4918345180071715, 1.2258868765909527, 0.9809228165973388, 0.6782876780174922, 0.23911709192210395, 0.0]),
array([3.885096429823169, 3.95183759232883, 4.306576404704602]),
array([3.3418364300367993, 1.6602456905723548, 1.2082676328391835, 0.7583238677565047, 0.5542220641369524]),
array([3.782815353284341, 3.140740668934252, 2.6232673646145064, 3.1929839760047165, 2.7434023612022846, 3.407688582135764, 2.9140902920144707, 2.6034529064722776, 3.4960690371096046, 3.5915136589282644, 3.1577455515147337, 3.5055220072876057, 2.836647963330244, 2.66716036287447, 1.2848865731627286, 1.6531837170525496, 1.1641455432459669, 1.5211604141998736, 1.3028968294036511, 1.716164501561674, 1.331683845699318, 1.487713726064366, 0.9619935458021238, 1.7419386771896819, 1.3734391868582219, 1.3328583788100645, 0.6708266054729306, 0.20503870873059948, 0.41613846380485203, 0.6294456831731645, 0.45191012690602345, 0.30187643439445916, 0.15343270415777088, 0.6366078672192506, 0.2005509392552789, 0.0]),
array([3.735916233089663, 3.157298607837066, 3.3361141250014503, 3.156908614613258, 2.862083225792999, 3.5907647551880215, 1.1112442275158372, 1.2544268036642996, 1.166148279665528, 1.0618031786308317, 0.9809071750247207, 1.501818278113653, 1.2888784299042921]),
array([2.998756930725609, 2.9744424814148736, 2.949557988166146, 3.353137478419505, 2.8716662221543108, 3.5067806839415883, 3.396744256803863, 2.7049805524124286, 2.8465023862589023, 3.585979534908428, 3.346656214422804, 3.381305387346972, 0.9439298680704553, 0.7917492975466074, 0.8935798220464053, 1.42147485938576, 1.6011253189541952, 1.301309004382686, 1.0442673605790729, 1.2980721169023899, 1.3076387002833718, 0.9458432009176452, 1.431729549130537, 0.923869226366933, 1.3854659078301936, 1.4655304653137722]),
array([3.039392013336707, 3.0540373447236715, 3.1706084336820073, 2.920867474544003, 2.748377910884716, 2.827759013428015, 2.7860908487832097, 3.1588760175221693, 2.564614805612462]),
array([2.6018291804204456, 1.3373842804875067, 1.023066549110309, 1.6716039336065898, 1.5998488190264144, 0.8304777226493092, 0.1690308052202939, 0.1807445433314966, 0.0]),
array([2.6219921232008137, 1.621174031698397, 1.2214858630282555, 1.2096785584714853, 0.36480359625544906, 0.3684674713452665, 0.5168970352645961, 0.0]),
array([0.2881793870910189, 0.0]),
array([1.3010726731708004, 1.1430045917002218, 1.565806566412227, 1.2368371491354937, 1.4207554529034119, 1.1257610624378565, 1.6560932093027607, 1.2182759122405735, 1.611209507864393, 1.5160627186857778, 1.5771347134702687, 1.6136883023136754, 0.4330024179801396, 0.0]),
array([0.9586005403596802, 1.5088639885535526, 0.8552608662801261, 0.8116448706732504, 1.2178855811340978, 1.6542818560036248, 1.4183353399962066, 0.3372780633454671, 0.500360072979039, 0.3018698160517439, 0.053928519336267314, 0.0]),
array([0.9871108236481698, 0.2409711559690204, 0.6850842519251554, 0.0]),
array([1.4256288474576853, 1.7113294545494584, 1.2326361757641764, 0.9388218523255405, 0.988951011500771, 0.23543103629888062, 0.0]),
array([1.7614729091713224, 1.6519917848292573, 1.7514890620441577]),
array([1.4838314264375971, 0.8942094757986385, 1.455547611486533, 1.443340635480308, 1.6794126017055917, 0.8972685383211596, 0.9674572629249879, 1.1160312873534806, 0.42710577167734876, 0.2277438295459031, 0.024515526508201727, 0.0]),
array([0.9905270240553113, 1.4924312519308516, 1.2354665041410242, 1.325326934014373, 1.5413214944703633, 0.49849738234663543, 0.532728482074112, 0.18163522020121936, 0.6164841297557971, 0.22804092693757627, 0.46396652500634034, 0.4442674851289249, 0.578953275624388, 0.40565735859454477, 0.5679866664390816, 0.0]),
array([1.2121365229046788, 0.0]),
array([0.8387807949002188, 0.7552785336717177, 0.6028409134994468, 0.0]),
array([1.4519986757835646, 1.2762688214678382, 1.3349288855922132, 1.3135354890876783, 0.8450042482363167, 1.1528197310902313, 1.1540567103550667, 0.9483576259237577, 1.3720099494005569, 1.1394015174840644, 1.1794403353293792, 1.0627594273078873, 0.7813091654279435, 0.9421552230366491, 0.7791436983686821, 0.4962987569664547, 0.49996756194245345, 0.6660098084355592]),
array([1.1476997868623868, 0.8860201911651054, 1.2860386067426994, 1.3714476150378432, 1.2654909411692576, 1.1826893623264159, 1.079056318707128, 1.0153863534104093, 0.8354411182164172, 0.3039170904075319, 0.5544970055411717, 0.0]),
array([0.8981448917174406, 0.3414339821086857, 0.2534433218048837, 0.0]),
array([0.16087465972658233, 0.5188609593184649, 0.2679398211637888, 0.6068168648359622, 0.030219722846439964, 0.0]),
array([0.13135590369831884, 0.6657236422845495, 0.0]),
array([0.27428022636858074, 0.10660861200858451, 0.0]),
array([0.25334476286240365, 0.24267096165427582, 0.4549750455328807, 0.0]),
array([0.14530164798943707, 0.0]),
array([0.3749202831687748, 0.24068296541508832, 0.0]),
array([0.354533232737813, 0.0])
]
d = [data_1]
names = ["42"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T2', 'T3', 'T4', 'T5', 'T6', 'T8', 'T9', 'T11', 'T12', 'T13', 'T14', 'T16', 'T17', 'T19', 'T20', 'T22', 'T23', 'T24', 'T25', 'T28', 'T29', 'T31', 'T32', 'T34', 'T35', 'T36', 'T37', 'T39', 'T41', 'T42', 'T43', 'T45', 'T46', 'T48', 'T49', 'T50', 'T52', 'T53', 'T54', 'T55', 'T56', 'T57', 'T58', 'T59', 'T61', 'T62', 'T63', 'T65', 'T66', 'T67', 'T69', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T80', 'T83', 'T86', 'T87', 'T90', 'T91', 'T92', 'T94', 'T95', 'T96', 'T98', 'T99', 'T102', 'T103', 'T104', 'T105', 'T108', 'T109', 'T110', 'T112', 'T114', 'T115', 'T116', 'T117', 'T118', 'T119', 'T120', 'T126', 'T128', 'T129', 'T131', 'T132', 'T134', 'T135', 'T136', 'T137', 'T138', 'T139', 'T142', 'T144', 'T146', 'T147', 'T148', 'T149', 'T151', 'T152', 'T153', 'T154', 'T156', 'T157', 'T158', 'T159', 'T160', 'T162', 'T163', 'T164', 'T166', 'T167', 'T169', 'T170', 'T171', 'T172', 'T173', 'T174', 'T176', 'T177', 'T178', 'T181', 'T182', 'T183', 'T184', 'T185', 'T186', 'T187', 'T188', 'T189', 'T191', 'T194', 'T195', 'T196', 'T197', 'T198', 'T200', 'T203', 'T204', 'T206', 'T207', 'T208', 'T209', 'T210']
def get_taxa_names(): return taxa_names