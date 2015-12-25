/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016 Klaus Spanderen

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_test_unix_parallel_test_runner_hpp
#define quantlib_test_unix_parallel_test_runner_hpp

#include <ql/types.hpp>

#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/thread/thread.hpp>
#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/interprocess/sync/named_mutex.hpp>

#define BOOST_TEST_NO_MAIN 1
#include <boost/test/results_collector.hpp>

#include <boost/test/included/unit_test.hpp>
#include <boost/algorithm/string.hpp>

#include <map>
#include <list>
#include <sstream>
#include <utility>
#include <fstream>
#include <iostream>

#include <string>
#include <cstring>
#include <cstdlib>
#include <unistd.h>

#ifdef BOOST_MSVC
#  error parallel test suite runner is not available on Windows
#endif

using boost::unit_test::test_results;
using namespace boost::interprocess;
using namespace boost::unit_test_framework;


namespace {
    class TestCaseCollector : public test_tree_visitor {
      public:
        typedef std::map<test_unit_id, std::list<test_unit_id> > id_map_t;

        const id_map_t& map() const { return idMap_; }
        test_unit_id testSuiteId() const { return testSuiteId_; }

        bool visit(test_unit const& tu) {
            if (tu.p_parent_id == framework::master_test_suite().p_id) {
                QL_REQUIRE(!tu.p_name.get().compare("QuantLib test suite"),
                     "could not find QuantLib test suite");
                testSuiteId_ = tu.p_id;
            }
            return test_tree_visitor::visit(tu);
        }

        void visit(test_case const& tc) {
            idMap_[tc.p_parent_id].push_back(tc.p_id);
        }

        std::list<test_unit_id>::size_type numberOfTests() {
            std::list<test_unit_id>::size_type n=0;
            for (id_map_t::const_iterator p_it = idMap_.begin();
                p_it != idMap_.end(); ++p_it) n+=p_it->second.size();

            return n;
        }
      private:
        id_map_t idMap_;
        test_unit_id testSuiteId_;
    };

    struct TestCaseId {
        test_unit_id id;
        bool terminate;
    };

    struct RuntimeLog {
        QuantLib::Time time;
        char testCaseName[256];
    };

    const char* const namesLogMutexName = "named_log_mutex";

    void output_logstream(
        std::ostream& out, std::streambuf* outBuf, std::stringstream& s) {

        static named_mutex mutex(open_or_create, "namesLogMutexName");
        scoped_lock<named_mutex> lock(mutex);

        out.flush();
        out.rdbuf(outBuf);

        std::vector<std::string> tok;
        const std::string lines = s.str();
        boost::split(tok, lines, boost::is_any_of("\n"));

        for (std::vector<std::string>::const_iterator iter = tok.begin();
            iter != tok.end(); ++iter) {
            if (iter->length() && iter->compare("Running 1 test case...")) {
                out << *iter << std::endl;
            }
        }

        s.str(std::string());
        out.rdbuf(s.rdbuf());
    }
}

test_suite* init_unit_test_suite(int, char* []);

int main( int argc, char* argv[] )
{
    typedef QuantLib::Time Time;

    const char* const profileFileName = ".unit_test_profile.txt";
    const char* const testUnitIdQueueName = "test_unit_queue";
    const char* const testResultQueueName = "test_result_queue";
    const char* const testRuntimeLogName  = "test_runtime_log_queue";

    std::map<std::string, Time> runTimeLog;

    std::ifstream in(profileFileName);
    if (in.good()) {
        for (std::string line; std::getline(in, line);) {
            std::vector<std::string> tok;
            boost::split(tok, line, boost::is_any_of(" "));

            QL_REQUIRE(tok.size() == 2,
                "every line should consists of two entries");
            runTimeLog[tok[0]] = std::atof(tok[1].c_str());
        }
    }
    in.close();

    unsigned int priority;
    message_queue::size_type recvd_size;

    try {
        framework::init(init_unit_test_suite, argc, argv );
        framework::finalize_setup_phase();

        TestCaseCollector tcc;
        traverse_test_tree(framework::master_test_suite(), tcc , true);

        s_log_impl().stream() << "Total number of test cases: "
            << tcc.numberOfTests() << std::endl;

        const unsigned nProc = boost::thread::hardware_concurrency();

        message_queue::remove(testUnitIdQueueName);
        message_queue mq(create_only, testUnitIdQueueName,
            tcc.numberOfTests() + nProc, sizeof(TestCaseId));

        message_queue::remove(testResultQueueName);
        message_queue rq(create_only, testResultQueueName, nProc,
            sizeof(test_results));

        message_queue::remove(testRuntimeLogName);
        message_queue lq(create_only, testRuntimeLogName, nProc,
            sizeof(RuntimeLog));

        // run root test cases in master process
        const std::list<test_unit_id>& qlRoot
            = tcc.map().find(tcc.testSuiteId())->second;

        test_results results;

        std::stringstream logBuf;
        std::streambuf* const oldBuf = s_log_impl().stream().rdbuf();
        s_log_impl().stream().rdbuf(logBuf.rdbuf());

        for (std::list<test_unit_id>::const_iterator iter = qlRoot.begin();
            std::distance(qlRoot.begin(), iter) < int(qlRoot.size())-1;
            ++iter) {
            framework::run(*iter);
            results += boost::unit_test::results_collector.results(*iter);
        }
        output_logstream(s_log_impl().stream(), oldBuf, logBuf);
        s_log_impl().stream().rdbuf(oldBuf);

        // fork worker processes
        pid_t pid;
        for (unsigned i = 0; (i < nProc) && (pid = fork()); ++i);

        if (pid) {
            struct mutex_remove {
                ~mutex_remove() { named_mutex::remove(namesLogMutexName); }
            } mutex_remover;

            struct queue_remove {
                queue_remove(const char* name) : name_(name) { }
                ~queue_remove() { message_queue::remove(name_); }

            private:
                const char* const name_;
            } queue_remover1(testUnitIdQueueName),
              queue_remover2(testResultQueueName),
              queue_remover3(testRuntimeLogName);

            std::multimap<Time, test_unit_id> testsSortedByRunTime;

            for (TestCaseCollector::id_map_t::const_iterator
                p_it = tcc.map().begin();
                p_it != tcc.map().end(); ++p_it) {

                if (p_it->first != tcc.testSuiteId()) {
                    for (std::list<test_unit_id>::const_iterator
                        it =  p_it->second.begin();
                        it != p_it->second.end(); ++it) {

                        const std::string& name
                            = framework::get(*it, TUT_ANY).p_name;

                        if (runTimeLog.count(name)) {
                            testsSortedByRunTime.insert(
                                std::make_pair(runTimeLog[name], *it));
                        }
                        else {
                            testsSortedByRunTime.insert(
                                std::make_pair(
                                    std::numeric_limits<Time>::max(), *it));
                        }
                    }
                }
            }

            std::list<test_unit_id> ids;
            for (std::multimap<Time, test_unit_id>::const_iterator
                iter = testsSortedByRunTime.begin();
                iter != testsSortedByRunTime.end(); ++iter) {
                ids.push_front(iter->second);
            }
            QL_REQUIRE(ids.size() + qlRoot.size() == tcc.numberOfTests(),
                "missing test case in distrubtion list");

            testsSortedByRunTime.clear();

            message_queue mq(open_only, testUnitIdQueueName);
            for (std::list<test_unit_id>::const_iterator iter = ids.begin();
                iter != ids.end(); ++iter) {
                const TestCaseId id = { *iter, false };
                mq.send(&id, sizeof(TestCaseId), 0);
            }

            const TestCaseId id = { 0, true };
            for (unsigned i = 0; i < nProc; ++i) {
                mq.send(&id, sizeof(TestCaseId), 0);
            }

            message_queue rq(open_only, testResultQueueName);
            for(unsigned i = 0; i < nProc; ++i) {
                test_results remoteResults;

                rq.receive(&remoteResults,
                    sizeof(remoteResults), recvd_size, priority);
                results+=remoteResults;
            }

            if (!qlRoot.empty()) {
                std::streambuf* const oldBuf = s_log_impl().stream().rdbuf();
                s_log_impl().stream().rdbuf(logBuf.rdbuf());

                framework::run(qlRoot.back());

                output_logstream(s_log_impl().stream(), oldBuf, logBuf);
                s_log_impl().stream().rdbuf(oldBuf);
            }

            results += boost::unit_test::results_collector.results(
                qlRoot.back());

            s_log_impl().stream() << "Test module \""
                << framework::master_test_suite().p_name
                <<"\" has passed with:"
                << std::endl
                << " " << results.p_assertions_failed << " assertions failed"
                << std::endl
                << " " << results.p_assertions_passed << " assertions passed"
                << std::endl;

            message_queue lq(open_only, testRuntimeLogName);

            RuntimeLog log;
            for (unsigned i=0; i < ids.size(); ++i) {
                lq.receive(&log, sizeof(RuntimeLog), recvd_size, priority);
                runTimeLog[std::string(log.testCaseName)] = log.time;
            }

            std::ofstream out(
                profileFileName, std::ios::out | std::ios::trunc);
            out << std::setprecision(6);
            for (std::map<std::string, QuantLib::Time>::const_iterator
                iter = runTimeLog.begin(); iter != runTimeLog.end(); ++iter) {
                out << iter->first << " " << iter->second << std::endl;
            }
            out.close();
        }
        else {
            std::stringstream logBuf;
            std::streambuf* const oldBuf = s_log_impl().stream().rdbuf();
            s_log_impl().stream().rdbuf(logBuf.rdbuf());

            message_queue mq(open_only, testUnitIdQueueName);

            TestCaseId id;
            mq.receive(&id, sizeof(TestCaseId), recvd_size, priority);

            typedef std::list<std::pair<std::string, QuantLib::Time> >
                run_time_list_type;
            run_time_list_type runTimeLogs;

            test_results results;

            while (!id.terminate) {

                boost::timer t;
                framework::run(id.id);

                runTimeLogs.push_back(std::make_pair(
                    framework::get(id.id, TUT_ANY).p_name, t.elapsed()));

                output_logstream(s_log_impl().stream(), oldBuf, logBuf);
                results+=boost::unit_test::results_collector.results(id.id);

                mq.receive(&id, sizeof(TestCaseId), recvd_size, priority);
            }

            message_queue rq(open_only, testResultQueueName);
            rq.send(&results, sizeof(results), 0);

            output_logstream(s_log_impl().stream(), oldBuf, logBuf);
            s_log_impl().stream().rdbuf(oldBuf);

            RuntimeLog log;
            log.testCaseName[sizeof(log.testCaseName)-1] = '\0';

            message_queue lq(open_only, testRuntimeLogName);
            for (run_time_list_type::const_iterator iter = runTimeLogs.begin();
                iter != runTimeLogs.end(); ++iter) {
                log.time = iter->second;

                std::strncpy(log.testCaseName, iter->first.c_str(),
                    sizeof(log.testCaseName)-1);

                lq.send(&log, sizeof(RuntimeLog), 0);
            }
        }
    }
    catch(QuantLib::Error &ex) {
        std::cerr << "QuantLib exception: " << ex.what() << std::endl;
        return boost::exit_exception_failure;
    }
    catch(interprocess_exception &ex){
        std::cerr << "interprocess exception: " << ex.what() << std::endl;
        return boost::exit_exception_failure;
    }
    catch( framework::nothing_to_test const& ) {
        return boost::exit_success;
    }
    catch( framework::internal_error const& ex ) {
        results_reporter::get_stream() << "Boost.Test framework internal error: " << ex.what() << std::endl;

        return boost::exit_exception_failure;
    }
    catch( framework::setup_error const& ex ) {
        results_reporter::get_stream() << "Test setup error: " << ex.what() << std::endl;

        return boost::exit_exception_failure;
    }
    catch( ... ) {
        results_reporter::get_stream() << "Boost.Test framework internal error: unknown reason" << std::endl;

        return boost::exit_exception_failure;
    }

    framework::shutdown();
}

#endif
