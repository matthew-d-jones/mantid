// WorkspaceFactory include must be first otherwise you get a bizarre Poco-related compilation error on Windows
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidICat/SearchHelper.h"
#include "MantidICat/Session.h"
#include "MantidICat/ErrorHandling.h" 
#include <iomanip>

namespace Mantid
{
	namespace ICat
	{
		using namespace Kernel;
		using namespace API;
		
		

		/* This method calls ICat API searchbydavanced and do the basic run search 
		 * @param icat Proxy object for ICat
		 * @param request request object
		 * @param response response object
		 */
		int CSearchHelper::doSearch(ICATPortBindingProxy& icat,boost::shared_ptr<ns1__searchByAdvanced>& request,ns1__searchByAdvancedResponse& response)
		{
			// Define ssl authentication scheme
			if (soap_ssl_client_context(&icat,
				SOAP_SSL_NO_AUTHENTICATION, /* use SOAP_SSL_DEFAULT in production code */
				NULL,       /* keyfile: required only when client must authenticate to
							server (see SSL docs on how to obtain this file) */
							NULL,       /* password to read the keyfile */
							NULL,      /* optional cacert file to store trusted certificates */
							NULL,      /* optional capath to directory with trusted certificates */
							NULL      /* if randfile!=NULL: use a file with random data to seed randomness */
							))
			{
				CErrorHandling::throwErrorMessages(icat);
			}
			clock_t start=clock();
			int ret_advsearch=icat.searchByAdvanced(request.get(),&response);
			if(ret_advsearch!=0)
			{
			
				CErrorHandling::throwErrorMessages(icat);
			}
			clock_t end=clock();
			float diff = float(end -start)/CLOCKS_PER_SEC;
			g_log.information()<<" Time taken to do  search is "<<diff<<std::endl;
			return ret_advsearch;
		}

		/* This method does a search  by run number and returns investigation data
		 * @param inputs reference to class containing search inputs
		 * @param responsews_sptr output table workspace
		 */
		int CSearchHelper::doIsisSearch(CSearchInput& inputs,API::ITableWorkspace_sptr& responsews_sptr)
		{
			//ICAt proxy object
			ICATPortBindingProxy icat;
			// request object
			boost::shared_ptr<ns1__searchByAdvanced> req_sptr(new ns1__searchByAdvanced );
			boost::shared_ptr<std::string > sessionId_sptr(new std::string);
			req_sptr->sessionId=sessionId_sptr.get();
			//get the sessionid which is cached in session class during login
			*req_sptr->sessionId=Session::Instance().getSessionId();

			boost::shared_ptr<ns1__advancedSearchDetails>adv_sptr(new ns1__advancedSearchDetails);
			req_sptr->advancedSearchDetails=adv_sptr.get();
			//run start
			boost::shared_ptr<double>runstart_sptr(new double);
			if(inputs.getRunStart()>0)
			{
				req_sptr->advancedSearchDetails->runStart=runstart_sptr.get();
			   *req_sptr->advancedSearchDetails->runStart=inputs.getRunStart();
			}
			//run end
			boost::shared_ptr<double>runend_sptr(new double);
			if(inputs.getRunEnd()>0)
			{
				req_sptr->advancedSearchDetails->runEnd=runend_sptr.get();
			    *req_sptr->advancedSearchDetails->runEnd=inputs.getRunEnd();
			}
            //start date
			boost::shared_ptr<time_t> startdate_sptr(new time_t);
			if(inputs.getStartDate()!=0)
			{				
				req_sptr->advancedSearchDetails->dateRangeStart = startdate_sptr.get();
				*req_sptr->advancedSearchDetails->dateRangeStart = inputs.getStartDate();
			}
			//end date
            boost::shared_ptr<time_t> enddate_sptr(new time_t);
			if(inputs.getEndDate()!=0)
			{				
				req_sptr->advancedSearchDetails->dateRangeEnd = enddate_sptr.get();
				*req_sptr->advancedSearchDetails->dateRangeEnd =inputs.getEndDate();
			}
			req_sptr->advancedSearchDetails->caseSensitive=inputs.getCaseSensitive();
			// investigation include
            boost::shared_ptr<ns1__investigationInclude>invstInculde_sptr(new ns1__investigationInclude);
			req_sptr->advancedSearchDetails->investigationInclude=invstInculde_sptr.get();
			// investigation include
			*req_sptr->advancedSearchDetails->investigationInclude=inputs.getInvestigationInclude();

			if(!inputs.getKeywords().empty())
			{
				req_sptr->advancedSearchDetails->keywords.push_back(inputs.getKeywords());
			}

			//setting the input parameters
			//setReqParamforSearchByRunNumber(inputs,req_sptr);
		
			//response object
			ns1__searchByAdvancedResponse response;
			// do  search
			int ret_search=doSearch(icat,req_sptr,response);
			if(ret_search!=0)
			{
				//replace with mantid error routine
				CErrorHandling::throwErrorMessages(icat);
			}
			if(response.return_.empty())
			{	
				g_log.information()<<"ICat investigations search is complete.There are no results to display"<<std::endl;
				return -1;
        		//throw std::runtime_error("ICat investigations search is complete.There are no results to display");
			}
			//save response to a table workspace
			saveSearchResponse(response,responsews_sptr);
			return ret_search;
		}
	   /**This method sets the input request parameters for search 
		* @param input Refrence to input class
		* @request refrence to request object
		*/
		void CSearchHelper::setReqParamforSearchByRunNumber(CSearchInput& input,boost::shared_ptr<ns1__searchByAdvanced>& request)
		{
			//get the sessionid which is cached in session class during login
			*request->sessionId=Session::Instance().getSessionId();
			//run start
			if(input.getRunStart()>0)
			{
			*request->advancedSearchDetails->runStart=input.getRunStart();
			}
			//run end
			if(input.getRunEnd()>0)
			{
			
			*request->advancedSearchDetails->runEnd=input.getRunEnd();
			}
			//instrument name
			if(!input.getInstrument().empty())
			{
				request->advancedSearchDetails->instruments.push_back(input.getInstrument());
			}
			if(!input.getKeywords().empty())
			{
				request->advancedSearchDetails->keywords.push_back(input.getKeywords());
			}

			if(input.getEndDate()!=0)
			{
			*request->advancedSearchDetails->dateRangeEnd =input.getEndDate();
			}
			if(input.getStartDate()!=0)
			{
				*request->advancedSearchDetails->dateRangeStart = input.getStartDate();
			}
			
			request->advancedSearchDetails->caseSensitive=input.getCaseSensitive();
			// investigation include
			*request->advancedSearchDetails->investigationInclude=input.getInvestigationInclude();
		}

	  		
	   /** This method saves the search response( investigations )data to a table workspace
		*  @param response const reference to response object
		*  @param outputws shared pointer to output workspace
		*  @returns shared pointer to table workspace which stores the data
		*/
		void  CSearchHelper::saveSearchResponse(const ns1__searchByAdvancedResponse& response,API::ITableWorkspace_sptr& outputws)
		{
			//create table workspace
		
			//API::ITableWorkspace_sptr outputws =createTableWorkspace();
			
			outputws->addColumn("long64","InvestigationId");
			outputws->addColumn("str","RbNumber");
			outputws->addColumn("str","Title");
			outputws->addColumn("str","Type");
			outputws->addColumn("str","Instrument");
			outputws->addColumn("str","Investigator");
			outputws->addColumn("str","RunRange");
			outputws->addColumn("str","Year");

			outputws->addColumn("str","Abstract");
			outputws->addColumn("str","Investigators First Name");
			outputws->addColumn("str","Investigators Second Name");
			outputws->addColumn("str","Samples Name");

			try
			{				
				saveInvestigations(response.return_,outputws);
			}
			catch(std::runtime_error& )
			{
			  throw std::runtime_error("Error when saving  the ICat Search Results data to Workspace");
			}
		
		}
	   /** This method saves investigations  to a table workspace
		*  @param investigations a vector containing investigation data
		*  @param outputws shared pointer to output workspace
		*/
		void CSearchHelper::saveInvestigations(const std::vector<ns1__investigation*>& investigations,API::ITableWorkspace_sptr& outputws)
		{
			
			try
			{
				std::vector<ns1__investigation*>::const_iterator citr;
				for (citr=investigations.begin();citr!=investigations.end();++citr)
				{
					API::TableRow t = outputws->appendRow();
					//investigation id
					savetoTableWorkspace((*citr)->id,t);
					//std::cout<<"investigation id is "<<(*(*citr)->id)<<std::endl;
										
					//rb number
					savetoTableWorkspace((*citr)->invNumber,t);

					//std::cout<<"rb number is "<<*(*citr)->invNumber<<std::endl;
					//title
					savetoTableWorkspace((*citr)->title,t);

					//std::cout<<"title  is "<<*(*citr)->title<<std::endl;
                   				
					//type 
					savetoTableWorkspace((*citr)->invType,t);

					savetoTableWorkspace((*citr)->instrument,t);
					//std::cout<<"instrument is "<<*(*citr)->instrument<<std::endl;
					//investigator
					savetoTableWorkspace((*citr)->bcatInvStr,t);
					// run range
					savetoTableWorkspace((*citr)->invParamValue,t);
								
					//year
					std::string *sInvEndtime=new std::string ;
					if((*citr)->invEndDate!=NULL)
					{
						time_t  invEndtime=*(*citr)->invEndDate;
						char temp [25];
						strftime (temp,25,"%H:%M:%S %Y-%d-%b",localtime(&invEndtime));
						std::string ftime(temp);
						
						sInvEndtime->assign(ftime);
						savetoTableWorkspace(sInvEndtime,t);
					}
					else
					{
						savetoTableWorkspace(sInvEndtime,t);//this is to write empty value to table workspace.
					}
                    // abstract
					savetoTableWorkspace((*citr)->invAbstract,t);

					std::vector<ns1__investigator*>investigators;
					investigators.assign((*citr)->investigatorCollection.begin(),(*citr)->investigatorCollection.end());
									
					std::string fullname; 
					//for loop for getting invetigator's first and last name
					std::vector<ns1__investigator*>::const_iterator invstrItr;
					for(invstrItr=investigators.begin();invstrItr!=investigators.end();++invstrItr)
					{
						std::string firstname;std::string lastname;std::string name;
						if((*invstrItr)->facilityUser)
						{
						
							if((*invstrItr)->facilityUser->firstName)
							{
								firstname = *(*invstrItr)->facilityUser->firstName;
							}
							if((*invstrItr)->facilityUser->lastName)
							{
								lastname = *(*invstrItr)->facilityUser->lastName;
							}
							name = firstname+" "+ lastname;
						}
						if(!fullname.empty())
						{
							fullname+=",";
						}
						fullname+=name;
					}//end of for loop for investigator's name.

					std::string* facilityUser = new std::string;
					facilityUser->assign(fullname);
                	//invetigator name
					savetoTableWorkspace(facilityUser,t);

					std::vector<ns1__sample*>samples;
					samples.assign((*citr)->sampleCollection.begin(),(*citr)->sampleCollection.end());
					std::string sNames;
					//for loop for samples name.
					std::vector<ns1__sample*>::const_iterator sItr;
					for(sItr=samples.begin();sItr!=samples.end();++sItr)
					{
						std::string sName;
						if((*sItr)->name)
						{
							//savetoTableWorkspace((*sItr)->name,t);
							sName=*((*sItr)->name);
						}
						if(!sNames.empty())
						{
							sNames+=",";
						}
						sNames+=sName;
					}
					std::string *samplenames = new std::string;
					samplenames->assign(sNames);
					savetoTableWorkspace(samplenames,t);
				}
			}
			catch(std::runtime_error& )
			{
			  throw std::runtime_error("Error when saving  the ICat Search Results data to Workspace");
			}
		}

		/* This method loops through the response return_vector and saves the datafile details to a table workspace
		 * @param response const reference to response object
		 * @returns shared pointer to table workspace which stores the data
		 */
		API::ITableWorkspace_sptr CSearchHelper::saveFileSearchResponse(const ns1__searchByAdvancedResponse& response)
		{
			//create table workspace
			API::ITableWorkspace_sptr outputws =createTableWorkspace();
			//add columns
			outputws->addColumn("str","Name");
			outputws->addColumn("int","File Size");
			outputws->addColumn("long64","FileId");
			outputws->addColumn("str","Format");
			outputws->addColumn("str","Format Version");
			outputws->addColumn("str","Format Type");
			outputws->addColumn("str","Create Time");

			std::vector<ns1__investigation*> investVec;
			investVec.assign(response.return_.begin(),response.return_.end());

			try
			{
				std::vector<ns1__investigation*>::const_iterator inv_citr;
				for (inv_citr=investVec.begin();inv_citr!=investVec.end();++inv_citr)
				{
					std::vector<ns1__dataset*> datasetVec;
					datasetVec.assign((*inv_citr)->datasetCollection.begin(),(*inv_citr)->datasetCollection.end());

					std::vector<ns1__dataset*>::const_iterator dataset_citr;
					for(dataset_citr=datasetVec.begin();dataset_citr!=datasetVec.end();++dataset_citr)
					{
						std::vector<ns1__datafile * >datafileVec;
						datafileVec.assign((*dataset_citr)->datafileCollection.begin(),(*dataset_citr)->datafileCollection.end());

						std::vector<ns1__datafile * >::const_iterator datafile_citr;
						for(datafile_citr=datafileVec.begin();datafile_citr!=datafileVec.end();++datafile_citr)
						{

							API::TableRow t = outputws->appendRow();
							savetoTableWorkspace((*datafile_citr)->name,t);
							savetoTableWorkspace((*datafile_citr)->fileSize,t);

							//long long fileId=*(*datafile_citr)->id;
							savetoTableWorkspace((*datafile_citr)->id,t);
							ns1__datafileFormat* fileFormat=(*datafile_citr)->datafileFormat;
							if(fileFormat)
							{
								if(fileFormat->datafileFormatPK)
								{
									savetoTableWorkspace((fileFormat->datafileFormatPK->name),t);
									savetoTableWorkspace((fileFormat->datafileFormatPK->version),t);
								}
								savetoTableWorkspace((fileFormat->formatType),t);

							}
							if((*datafile_citr)->datafileCreateTime!=NULL)
							{
								time_t  crtime=*(*datafile_citr)->datafileCreateTime;
								char temp [25];
								strftime (temp,25,"%H:%M:%S %Y-%d-%b",localtime(&crtime));
								std::string ftime(temp);
								std::string *creationtime=new std::string ;
								creationtime->assign(ftime);
								savetoTableWorkspace(creationtime,t);
							}
			
						}

					}


				}
			}

			catch(std::runtime_error& )
			{
				throw;
			}

			return outputws;
		}
		/* This method sets the request parameters for the investigations includes
		 * @param invstId - investigation id 
		 * @param include - enum parameter to retrieve dat from DB
		 * @param request - request object
		*/
		void CSearchHelper::setReqParamforInvestigationIncludes(long long invstId,ns1__investigationInclude include,ns1__getInvestigationIncludes& request)
		{
			//get the sessionid which is cached in session class during login
			*request.sessionId=Session::Instance().getSessionId();;
			*request.investigationInclude=include;
  		    *request.investigationId=invstId;

		}
		/**This method calls ICat API getInvestigationIncludes and returns investigation details for a given investigation Id
		  *@param invstId - investigation id
		  *@param include - enum parameter for selecting the response data from the db.
		  *@param responsews_sptr - table workspace to save the response data
		*/
		int CSearchHelper::getDataFiles(long long invstId,bool bDataFiles,ns1__investigationInclude include,
			               API::ITableWorkspace_sptr& responsews_sptr)
		{
			//ICAt proxy object
			ICATPortBindingProxy icat;
			// Define ssl authentication scheme
			if (soap_ssl_client_context(&icat,
				SOAP_SSL_NO_AUTHENTICATION, /* use SOAP_SSL_DEFAULT in production code */
				NULL,       /* keyfile: required only when client must authenticate to
							server (see SSL docs on how to obtain this file) */
							NULL,       /* password to read the keyfile */
							NULL,      /* optional cacert file to store trusted certificates */
							NULL,      /* optional capath to directory with trusted certificates */
							NULL      /* if randfile!=NULL: use a file with random data to seed randomness */
							))
			{
				CErrorHandling::throwErrorMessages(icat);
			}

			ns1__getInvestigationIncludes request;
			//get the sessionid which is cached in session class during login
			boost::shared_ptr<std::string >sessionId_sptr(new std::string);
			request.sessionId=sessionId_sptr.get();

			// enum include
			boost::shared_ptr<ns1__investigationInclude>invstInculde_sptr(new ns1__investigationInclude);
			request.investigationInclude=invstInculde_sptr.get();

			//boost::shared_ptr<long long>invstId_sptr(new long long);
		   // request.investigationId=invstId_sptr.get();
			request.investigationId=new LONG64;
			setReqParamforInvestigationIncludes(invstId,include,request);

			ns1__getInvestigationIncludesResponse response;
			int ret_advsearch=icat.getInvestigationIncludes(&request,&response);
			if(ret_advsearch!=0)
			{
				CErrorHandling::throwErrorMessages(icat);
			}
			std::stringstream stream;
			stream<<invstId;
			if(!response.return_)
			{
				//throw std::runtime_error("No data files exists in the ICat database for the selected investigation");
				g_log.information()<<"No data files exists in the ICat database for the selected investigation"<<std::endl;
				return -1;
			}
			try
			{
				//responsews_sptr=saveInvestigationIncludesResponse(bDataFiles,response);
				saveInvestigationIncludesResponse(bDataFiles,response,responsews_sptr);
			}
			catch(std::runtime_error)
			{				
				throw std::runtime_error("Error when selecting the investigation data with inestigation id "+ stream.str());
			}
			return ret_advsearch;
		}
		/* This method loops through the response return_vector and saves the datafile details to a table workspace
		 * @param response const reference to response object
		 * @param outputws shared pointer to table workspace which stores the data
		*/
		void  CSearchHelper::saveInvestigationIncludesResponse(bool bloadonlyData,
			const ns1__getInvestigationIncludesResponse& response,API::ITableWorkspace_sptr& outputws)
		{
			//create table workspace
			//API::ITableWorkspace_sptr outputws =createTableWorkspace();

			//outputws->addColumn("str","Instrument");//Instrument name
			//outputws->addColumn("long64","InvestigationId");//investigation id
			outputws->addColumn("str","Name");//File name
			outputws->addColumn("int","File Size");//File Size
			outputws->addColumn("long64","File Id");//File id
			outputws->addColumn("str","Format");//File Format
			outputws->addColumn("str","Format Version");//File Version
			outputws->addColumn("str","Format Type");// File Format Type
			outputws->addColumn("str","Create Time");// File Creation Time

			try
			{		std::vector<ns1__dataset*> datasetVec;
					datasetVec.assign((response.return_)->datasetCollection.begin(),(response.return_)->datasetCollection.end());
					if(datasetVec.empty())
					{
						throw std::runtime_error("No data files exists in the database for the selected investigation");
					}
					std::vector<ns1__dataset*>::const_iterator dataset_citr;
					for(dataset_citr=datasetVec.begin();dataset_citr!=datasetVec.end();++dataset_citr)
					{
						std::vector<ns1__datafile * >datafileVec;
						datafileVec.assign((*dataset_citr)->datafileCollection.begin(),(*dataset_citr)->datafileCollection.end());
						if(datafileVec.empty())
						{
							throw std::runtime_error("No data files exists in the database for the selected  investigation ");
						}

						std::vector<ns1__datafile * >::const_iterator datafile_citr;
						for(datafile_citr=datafileVec.begin();datafile_citr!=datafileVec.end();++datafile_citr)
						{

                           if(bloadonlyData)
						   {
							if(!isDataFile((*datafile_citr)->name))
							{
								continue;
							}
						   }
							API::TableRow t = outputws->appendRow();
							
							//instrument name
							//savetoTableWorkspace(response.return_->instrument,t);
							//investigation Id
							//savetoTableWorkspace(response.return_->id,t);
							// File Name
							savetoTableWorkspace((*datafile_citr)->name,t);
						    // File Size
							savetoTableWorkspace((*datafile_citr)->fileSize,t);
							//File Id
							savetoTableWorkspace((*datafile_citr)->id,t);
							ns1__datafileFormat* fileFormat=(*datafile_citr)->datafileFormat;
							if(fileFormat)
							{
								if(fileFormat->datafileFormatPK)
								{
									// File Format
									savetoTableWorkspace((fileFormat->datafileFormatPK->name),t);
									// File Format Version
									savetoTableWorkspace((fileFormat->datafileFormatPK->version),t);
								}
								else
								{
									std::string *s=NULL;
									savetoTableWorkspace(s,t);
									savetoTableWorkspace(s,t);
								}
								// File format Type
								savetoTableWorkspace((fileFormat->formatType),t);
							}
							else
							{
								//i've to see a better way to write empty columns if the data is empty
								std::string *s=NULL;
								savetoTableWorkspace(s,t);
								savetoTableWorkspace(s,t);
								savetoTableWorkspace(s,t);
							}
														
							//File creation Time.
							std::string *creationtime=NULL;
							if((*datafile_citr)->datafileCreateTime!=NULL)
							{
								time_t  crtime=*(*datafile_citr)->datafileCreateTime;
								char temp [25];
								strftime (temp,25,"%H:%M:%S %Y-%d-%b",localtime(&crtime));
								std::string ftime(temp);
								creationtime=new std::string ;
								creationtime->assign(ftime);
							}
							savetoTableWorkspace(creationtime,t);

						}

					}

			}
			catch(std::runtime_error& )
			{
				throw ;
			}

			//return outputws;
		}

		/**This checks the datafile boolean  selected
		 * @param fileName - pointer to file name
		 * @return bool - returns true if it's a raw file or nexus file
		*/

		bool CSearchHelper::isDataFile(const std::string* fileName)
		{	
			if(!fileName)
			{
				return false;
			}
			std::basic_string <char>::size_type dotIndex;
			const std::basic_string <char>::size_type npos = -1;
			//find the position of .in row file
			dotIndex = (*fileName).find_last_of (".");
			std::string fextn=(*fileName).substr(dotIndex+1,(*fileName).size()-dotIndex);
			bool bData;
			(!fextn.compare("raw")|| !fextn.compare("RAW")|| !fextn.compare("nxs")|| !fextn.compare("NXS")) ? bData = true : bData = false;
			return bData;
		}

		/**This method calls ICat API getInvestigationIncludes and returns datasets details for a given investigation Id
		 * @param invstId - investigation id
		 * @param include - enum parameter for selecting the response data from iact db.
		 * @param responsews_sptr - table workspace to save the response data
		*/
		int CSearchHelper::doDataSetsSearch(long long invstId,ns1__investigationInclude include,
			               API::ITableWorkspace_sptr& responsews_sptr)
		{
			//ICAt proxy object
			ICATPortBindingProxy icat;
			// Define ssl authentication scheme
			if (soap_ssl_client_context(&icat,
				SOAP_SSL_NO_AUTHENTICATION, /* use SOAP_SSL_DEFAULT in production code */
				NULL,       /* keyfile: required only when client must authenticate to 
							server (see SSL docs on how to obtain this file) */
							NULL,       /* password to read the keyfile */
							NULL,      /* optional cacert file to store trusted certificates */
							NULL,      /* optional capath to directory with trusted certificates */
							NULL      /* if randfile!=NULL: use a file with random data to seed randomness */ 
							))
			{ 
				CErrorHandling::throwErrorMessages(icat);
			}

			// request object
			ns1__getInvestigationIncludes request;
			//get the sessionid which is cached in session class during login
			boost::shared_ptr<std::string >sessionId_sptr(new std::string);
			request.sessionId=sessionId_sptr.get();

			// enum include 
			boost::shared_ptr<ns1__investigationInclude>invstInculde_sptr(new ns1__investigationInclude);
			request.investigationInclude=invstInculde_sptr.get();

			request.investigationId=new LONG64;
			setReqParamforInvestigationIncludes(invstId,include,request);
            
			//response object
			ns1__getInvestigationIncludesResponse response;
			// Calling Icat api 
			int ret_advsearch=icat.getInvestigationIncludes(&request,&response);
			if(ret_advsearch!=0)
			{
				CErrorHandling::throwErrorMessages(icat);
			}
			std::stringstream stream;
			stream<<invstId;
			if(!(response.return_)|| (response.return_)->datasetCollection.empty())
			{
				//throw std::runtime_error("No datasets  exists in the ICat database for the inestigation id "+ stream.str());
				g_log.information()<<"No datasets  exists in the ICat database for the inestigation id "+ stream.str()<<std::endl;
				return -1 ;
			}
			try
			{
				//responsews_sptr=saveDataSets(response);
				saveDataSets(response,responsews_sptr);
			}
			catch(std::runtime_error)
			{
				
				throw std::runtime_error("Error when loading the datasets for the inestigation id "+ stream.str());
			}

			return ret_advsearch;
		}

		/* This method loops through the response return_vector and saves the datasets details to a table workspace
		 * @param response const reference to response object
		 * @returns shared pointer to table workspace which stores the data
		*/
		void  CSearchHelper::saveDataSets(const ns1__getInvestigationIncludesResponse& response,API::ITableWorkspace_sptr& outputws)
		{
			//create table workspace
			//API::ITableWorkspace_sptr outputws =createTableWorkspace();
			//adding columns
			outputws->addColumn("str","Name");//File name
			outputws->addColumn("str","Status");
			outputws->addColumn("str","Type");
			outputws->addColumn("str","Description");
			outputws->addColumn("long64","Sample");
			
			try
			{	
				
				std::vector<ns1__dataset*> datasetVec;
				datasetVec.assign((response.return_)->datasetCollection.begin(),(response.return_)->datasetCollection.end());

				std::vector<ns1__dataset*>::const_iterator dataset_citr;
				for(dataset_citr=datasetVec.begin();dataset_citr!=datasetVec.end();++dataset_citr)
				{						
					API::TableRow t = outputws->appendRow();

					// DataSet Name
					savetoTableWorkspace((*dataset_citr)->name,t);
					// DataSet Status
					savetoTableWorkspace((*dataset_citr)->datasetStatus,t);
					//DataSet Type
					savetoTableWorkspace((*dataset_citr)->datasetType,t);
					//DataSet Type
					savetoTableWorkspace((*dataset_citr)->description,t);
					//DataSet Type
					savetoTableWorkspace((*dataset_citr)->sampleId,t);
				}
			}

			catch(std::runtime_error& )
			{
				throw;
			}

			//return outputws;
		}

		/* *This method calls ICat api listruments and returns the list of instruments a table workspace
		   *@return API::ITableWorkspace_sptr
		*/
		API::ITableWorkspace_sptr CSearchHelper::listInstruments()
		{		
			//ICAt proxy object
			ICATPortBindingProxy icat;
			// Define ssl authentication scheme
			if (soap_ssl_client_context(&icat,
				SOAP_SSL_NO_AUTHENTICATION, /* use SOAP_SSL_DEFAULT in production code */
				NULL,       /* keyfile: required only when client must authenticate to 
							server (see SSL docs on how to obtain this file) */
							NULL,       /* password to read the keyfile */
							NULL,      /* optional cacert file to store trusted certificates */
							NULL,      /* optional capath to directory with trusted certificates */
							NULL      /* if randfile!=NULL: use a file with random data to seed randomness */ 
							))
			{ 
				CErrorHandling::throwErrorMessages(icat);
			}

			ns1__listInstruments request;

			//get the sessionid which is cached in session class during login
			boost::shared_ptr<std::string >sessionId_sptr(new std::string);
			request.sessionId=sessionId_sptr.get();
			//setting the request parameters
			setReqparamforlistInstruments(request);

			// response object
			 ns1__listInstrumentsResponse response;

			 int ret=icat.listInstruments(&request,&response);
			 if(ret!=0)
			 {
				 CErrorHandling::throwErrorMessages(icat);
			 }
			 if(response.return_.empty())
			 {
				 std::runtime_error("There are no instruments entry in the ICat data base");
			 }
	
			return saveInstrumentList(response);
		}

		/**This method sets the request parameter for ICat api list isnturments
		 * @param request - reference to request object
		*/
		void CSearchHelper::setReqparamforlistInstruments(ns1__listInstruments& request)
		{
			*request.sessionId=Session::Instance().getSessionId();
		}

		/**This method saves the response data for ICat api list isnturments
		 * @param response - reference to response object
		 * @return API::ITableWorkspace_sptr - shared pointer to table workspace
		 */
		API::ITableWorkspace_sptr  CSearchHelper::saveInstrumentList(const ns1__listInstrumentsResponse& response)
		{			
			API::ITableWorkspace_sptr outputws =createTableWorkspace();
			outputws->addColumn("str","Instrument Name");
			try
			{			
				std::vector<std::string>::const_iterator inst_citr;
				for(inst_citr=response.return_.begin();inst_citr!=response.return_.end();++inst_citr)
				{						
					API::TableRow t = outputws->appendRow();
					// instrument  Name
					t<<*inst_citr;
				}
			}
			catch(std::runtime_error& )
			{
				throw;
			}

			return outputws;
		}

		/// This method creates table workspace
		API::ITableWorkspace_sptr CSearchHelper::createTableWorkspace()
		{
			//create table workspace
			API::ITableWorkspace_sptr outputws ;
			try
			{
				outputws=WorkspaceFactory::Instance().createTable("TableWorkspace");
			}
			catch(Mantid::Kernel::Exception::NotFoundError& )
			{
				throw std::runtime_error("Error when saving  the ICat Search Results data to Workspace");
			}
			catch(std::runtime_error)
			{
				throw std::runtime_error("Error when saving  the ICat Search Results data to Workspace");
			}
			return outputws;
		}
        /// This method calls ICat api logoutand disconnects from ICat DB
		int CSearchHelper::doLogout()
		{
			ICATPortBindingProxy icat;
			// Define ssl authentication scheme
			if (soap_ssl_client_context(&icat,
				SOAP_SSL_NO_AUTHENTICATION, /* use SOAP_SSL_DEFAULT in production code */
				NULL,       /* keyfile: required only when client must authenticate to
							server (see SSL docs on how to obtain this file) */
							NULL,       /* password to read the keyfile */
							NULL,      /* optional cacert file to store trusted certificates */
							NULL,      /* optional capath to directory with trusted certificates */
							NULL      /* if randfile!=NULL: use a file with random data to seed randomness */
							))
			{
				
				CErrorHandling::throwErrorMessages(icat);
			}

			ns1__logout request;
			ns1__logoutResponse response;
			boost::shared_ptr<std::string > sessionId_sptr(new std::string);
			*sessionId_sptr=Session::Instance().getSessionId();
			request.sessionId=sessionId_sptr.get();
			int ret=icat.logout(&request,&response);
			if(ret!=0)
			{
				//icat.soap_stream_fault(std::cerr);
				CErrorHandling::throwErrorMessages(icat);
			}
			
			return ret;
		}

		/**This method calls ICat api getmyinvestigations and do returns the investigations of the logged in user
		 * @param ws_sptr - shared pointer to table workspace which stores the investigations search result
		 */
		void CSearchHelper::doMyDataSearch(API::ITableWorkspace_sptr& ws_sptr)
		{
			ICATPortBindingProxy icat;
			// Define ssl authentication scheme
			if (soap_ssl_client_context(&icat,
				SOAP_SSL_NO_AUTHENTICATION, /* use SOAP_SSL_DEFAULT in production code */
				NULL,       /* keyfile: required only when client must authenticate to
							server (see SSL docs on how to obtain this file) */
							NULL,       /* password to read the keyfile */
							NULL,      /* optional cacert file to store trusted certificates */
							NULL,      /* optional capath to directory with trusted certificates */
							NULL      /* if randfile!=NULL: use a file with random data to seed randomness */
							))
			{
				
				CErrorHandling::throwErrorMessages(icat);
			}

			ns1__getMyInvestigationsIncludes request;
			ns1__getMyInvestigationsIncludesResponse response;
			boost::shared_ptr<std::string > sessionId_sptr(new std::string);
			*sessionId_sptr=Session::Instance().getSessionId();
			request.sessionId=sessionId_sptr.get();
			// investigation include
            boost::shared_ptr<ns1__investigationInclude>invstInculde_sptr(new ns1__investigationInclude);
			request.investigationInclude=invstInculde_sptr.get();
			*request.investigationInclude = ns1__investigationInclude__INVESTIGATORS_USCOREONLY;

			int ret=icat.getMyInvestigationsIncludes(&request,&response);
			if(ret!=0)
			{
				CErrorHandling::throwErrorMessages(icat);
			}
			if(response.return_.empty())
			{	
				g_log.information()<<"ICat Mydata search is complete.There are no results to display"<<std::endl;
				return ;
        		//throw std::runtime_error("ICat investigations search is complete.There are no results to display");
			}
			//save response to a table workspace
			saveMyInvestigations(response,ws_sptr);
			
			
		}

		/**This method calls ICat api getmyinvestigations and do returns the investigations of the logged in user
		 * @param response - reference to response  object 
		 * @param outputws - shared pointer to table workspace which stores the investigations search result
		 */
		void CSearchHelper::saveMyInvestigations( const ns1__getMyInvestigationsIncludesResponse& response,API::ITableWorkspace_sptr& outputws)
		{
			outputws->addColumn("long64","InvestigationId");
			outputws->addColumn("str","RbNumber");
			outputws->addColumn("str","Title");
			outputws->addColumn("str","Type");
			outputws->addColumn("str","Instrument");
			outputws->addColumn("str","Investigator");
			outputws->addColumn("str","RunRange");
			outputws->addColumn("str","Year");
			
			outputws->addColumn("str","Abstract");
			outputws->addColumn("str","Investigators First Name");
			outputws->addColumn("str","Investigators Second Name");
			outputws->addColumn("str","Samples Name");

			saveInvestigations(response.return_,outputws);

		}




	}
}
