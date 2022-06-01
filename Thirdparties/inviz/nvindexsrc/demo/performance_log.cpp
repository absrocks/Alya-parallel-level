/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief NVIDIA IndeX demo viewer performance logging functions

#include "performance_log.h"

#include "common/string_dict.h"
#include "common/tokenizer.h"

#include "nvindex_appdata.h"
#include "nvindex_rendering_context.h"

#include <iterator>

//----------------------------------------------------------------------
// lock for the logger set 
static mi::base::Lock Performance_logger_set_lock;
static mi::base::Lock Performance_logger_logging_set_state_lock;

//----------------------------------------------------------------------
void update_recent_n_statistics(
    nv::index::IPerformance_values * performance_values)
{
    if(Nvindex_AppData::instance()->is_enabled_recent_n_stat())
    {
        // cumulative fps statistics: max and average
        const mi::Float64 frames_per_sec = performance_values->get_time("frames_per_second");
        Cumulative_stat * p_cstat = Nvindex_AppData::instance()->peek_cumulative_stat();
        p_cstat->add_sample(frames_per_sec);
        Nvindex_AppData::instance()->set_stat_recent_n_max_fps(static_cast<mi::Float32>(p_cstat->get_max()));
        Nvindex_AppData::instance()->set_stat_recent_n_ave_fps(static_cast<mi::Float32>(p_cstat->get_ave()));
    }
}

//----------------------------------------------------------------------
void set_monitor_performance_values(mi::base::Handle<nv::index::IConfig_settings> & edit_config_settings,
                                    bool is_monitor_on)
{
    assert(edit_config_settings.is_valid_interface());
    edit_config_settings->set_monitor_performance_values(is_monitor_on);

    Performance_logger_app_state & perf_log_state = Nvindex_AppData::instance()->peek_log_state();
    perf_log_state.set_monitoring_on(is_monitor_on);
}

//----------------------------------------------------------------------
void set_performance_monitoring_by_project(
    Nvindex_rendering_context& irc,
    nv::index_common::String_dict & app_proj,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    assert(dice_transaction != 0);

    // get the parameters
    const bool is_monitor = 
        nv::index_common::get_bool(app_proj.get("index::config::set_monitor_performance_values", "no"));

    // performance_log_state and initialization lock (peek_log_state())
    mi::base::Lock::Block block(&Performance_logger_set_lock);

    // Initialize the performance log application state
    Performance_logger_app_state & perf_log_state = Nvindex_AppData::instance()->peek_log_state();
    for(mi::Uint32 i = 0; i < Performance_logger_base::LK_COUNT; ++i){
        const std::string key_name =  Performance_logger_base::get_log_type_string(i);
        const std::string enable_key = "app::performance::" + key_name + "::enable_log";
        const std::string file_key   = "app::performance::" + key_name + "::file";
        const bool is_logging        =  nv::index_common::get_bool(app_proj.get(enable_key, "no"));
        const std::string fname      =  app_proj.get(file_key, "/dev/null");
            
        perf_log_state.set_logging_on(i, is_logging);
        perf_log_state.set_logging_filename(i, fname);
        if (is_logging){
            INFO_LOG << "perf log: " << enable_key << ": " << is_logging << ", " << file_key << ": " << fname;
        }
    }
    if(!is_monitor){
        INFO_LOG << "Performance monitoring is currently switched off.";
    }

    // we got the all performance logging parameters now.
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(irc.m_session_tag));
    assert(session.is_valid_interface());
    mi::base::Handle<nv::index::IConfig_settings> config_settings(
        dice_transaction->edit<nv::index::IConfig_settings>(session->get_config()));
    set_monitor_performance_values(config_settings, is_monitor); // for sync state

    // Create the performance logger
    {
        Performance_logger_set & perf_logger_set = Nvindex_AppData::instance()->peek_logger_set();
        perf_logger_set.initialize_logger(app_proj.get("app::performance::logger_type", "csv"),
                                          app_proj);
    }
}


//----------------------------------------------------------------------
/// update performance monitoring entry
///
/// \param[in] p_perf_logger     performance logger
/// \param[in] performance_value an IPerformance_values
/// \param[in] app_proj          application project
/// \param[in] perf_logger_type  performance logger type (system, per_host, per_span)
static void update_performance_logging_entry(
    Performance_logger_base * p_perf_logger,
    nv::index::IPerformance_values * performance_value,
    nv::index_common::String_dict & app_proj,
    mi::Uint32 perf_logger_type)
{
    // check it has already been updated
    assert(p_perf_logger != 0);
    const std::vector<std::string> & logging_entry = p_perf_logger->peek_performance_logging_entry();
    if(!(logging_entry.empty())){
        // it has already been updated
        return;
    }

    // needs update
    const std::string key_name =  Performance_logger_base::get_log_type_string(perf_logger_type);
    const std::string item_list_key = "app::performance::" + key_name + "::item_list";
    const std::string item_list_str = app_proj.get(item_list_key, "<all>");
    std::vector<std::string> item_vec;
    if(item_list_str == "<all>"){
        INFO_LOG << "performance item_list is <all>, expand it.";
        if((perf_logger_type == Performance_logger_base::LK_System_global) ||
           (perf_logger_type == Performance_logger_base::LK_Per_host))
        {
            item_vec = Performance_logger_base::get_per_host_performance_value_entry_from_performance_value(
                performance_value);
        }
        else if(perf_logger_type == Performance_logger_base::LK_Per_span){
            item_vec = Performance_logger_base::get_per_span_performance_value_entry_from_performance_value(
                performance_value);
        }
        else{
            ERROR_LOG << "invalid perf_logger_type [" << perf_logger_type << "], use LK_System_global.";
            item_vec = Performance_logger_base::get_per_host_performance_value_entry_from_performance_value(
                performance_value);
        }
    }
    else{
        // just set the list from the project
        nv::index_common::Tokenizer::parse(item_list_str, " ", item_vec);
    }
    std::stringstream sstr;
    std::copy(item_vec.begin(), item_vec.end(), std::ostream_iterator<std::string>(sstr, ", "));
    INFO_LOG << "performance_logging_entry for [" << key_name << "]: " << sstr.str();

    p_perf_logger->set_performance_logging_entry(item_vec);
}

//----------------------------------------------------------------------
/// \return true p_logger_ref's file is a disk file and it exists.
static bool is_need_reopen(Performance_logger_base * p_logger_ref)
{
    assert(p_logger_ref != 0);
    
    if(!p_logger_ref->is_stream_file_on_disk()){
        // stdout: no need to reopen
        return false;
    }

    if(!Filename_ostream::isfile(p_logger_ref->get_filename())){
        // the file is a disk file, but not exist
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
void log_performance_information(
    nv::index::IPerformance_values * performance_value,
    mi::Sint32 frame_num)
{
    // check the logger need to be initialized 

    // performance_log_state and initialization lock (peek_log_state())
    mi::base::Lock::Block block(&Performance_logger_set_lock);

    // open close: compare the current logger state and log_app_state
    const Performance_logger_app_state & perf_log_state = Nvindex_AppData::instance()->peek_log_state();
    if(!(perf_log_state.is_monitoring_on())){
        return;                 // no performance monitoring
    }

    Performance_logger_set & logger_set = Nvindex_AppData::instance()->peek_logger_set();
    if(!(logger_set.is_valid())){
        ERROR_LOG << "invalid logger_set. No performance logging.";
        return;
    }

    nv::index_common::String_dict * p_app_proj = Nvindex_AppData::instance()->peek_app_proj();

    // reopen the stream if needed.
    for(mi::Uint32 i = 0; i < Performance_logger_base::LK_COUNT; ++i){
        if(perf_log_state.is_logging_on(i)){
            Performance_logger_base * p_logger_ref = logger_set.peek_logger(i);
            assert(p_logger_ref != 0);
            // if the stream is not good, need reopen.
            if((!(p_logger_ref->good())) || is_need_reopen(p_logger_ref)){

                // set the current filename
                p_logger_ref->set_filename(perf_log_state.get_logging_filename(i));
                if(is_need_reopen(p_logger_ref)){
                    INFO_LOG << "Lost the performance output file [" << p_logger_ref->get_filename()
                             << "], reopen it.";
                }

                // then open
                if((p_logger_ref->reopen_stream())){
                    // When re-open the stream, we write the header.
                    update_performance_logging_entry(p_logger_ref, performance_value, *p_app_proj, i);
                    p_logger_ref->append_header(*p_app_proj);
                }
                else{
                    ERROR_LOG << "Fail to open the log stream.";
                }
            }
        }
    }

    // log performance value when enabled
    for(mi::Uint32 i = 0; i < Performance_logger_base::LK_COUNT; ++i){
        if(perf_log_state.is_logging_on(i)){
            Performance_logger_base * p_logger_ref = logger_set.peek_logger(i);
            assert(p_logger_ref != 0);
            assert(p_logger_ref->good());
            p_logger_ref->append_performance_log(performance_value, frame_num);
        }
    }
}

//----------------------------------------------------------------------
bool setup_performance_logger_by_command_str(const std::string & command, std::string & ret_mes)
{
    std::string separators = " ";
    std::vector< std::string > tokens;
    nv::index_common::Tokenizer::parse(command, separators, tokens);
    
    const mi::Size token_size = tokens.size();
    if((token_size <= 1) || (tokens[0] != "perf_logger")){
        ret_mes = "Invalid command: perf_logger [logger_type|?] [option]";
        return false;
    }

    // Now, token_size >= 2, tokens[0] == "perf_log".
    Performance_logger_set & logger_set = Nvindex_AppData::instance()->peek_logger_set();
    if(tokens[1] == "?"){
        ret_mes = "current logger type [" + logger_set.get_logger_type_name() + "]";
        INFO_LOG << "perf_logger [logger_type] [logger_option]\n"
                 << "    logger_type: csv, human_readable\n"
                 << "    logger_option: simple, transfer (only effective when human_readable)\n"
                 << "Examples:\n"
                 << "     perf_logger ?\n"
                 << "     perf_logger csv\n"
                 << "     perf_logger human_readable simple\n"
                 << "     perf_logger human_readable transfer\n"
                 << "     perf_logger qa\n"
                 << ret_mes;
        return true;
    }
    else if((tokens[1] == "csv") || (tokens[1] == "qa")){
        const std::string logger_name = tokens[1];

        if(logger_set.get_logger_type_name() == logger_name){
            ret_mes = "Current perf_logger is " + logger_name + ", no change.";            
            return true;
        }
        nv::index_common::String_dict dummy_dict;
        {
            // Lock the logger set when delete and new them
            mi::base::Lock::Block block(&Performance_logger_set_lock);
            logger_set.initialize_logger(logger_name, dummy_dict);
        }
        ret_mes = "Set perf_logger to " + logger_name + ".";
        return true;
    }
    else if(tokens[1] == "human_readable"){
        nv::index_common::String_dict log_dict;
        std::string human_readable_type = "simple";
        if(token_size == 2){
            // use simple
        }
        else if(token_size == 3){
            if((tokens[2] == "simple") || (tokens[2] == "transfer")){
                human_readable_type = tokens[2];
            }
            else{
                ret_mes = "Unknown human_readable_type: " + tokens[2];
                return false;
            }
        }
        else{
            ret_mes = "Too many command args [" + command + "]";
            return false;
        }

        log_dict.insert("app::performance::global::human_readable_type",   human_readable_type);
        log_dict.insert("app::performance::per_host::human_readable_type", human_readable_type);
        log_dict.insert("app::performance::per_span::human_readable_type", human_readable_type);
        {
            // Lock the logger set when delete and new them
            mi::base::Lock::Block block(&Performance_logger_set_lock);
            logger_set.initialize_logger("human_readable", log_dict);
        }
        ret_mes = "Set perf_logger to human_readable " + tokens[2];
        return true;
    }

    ret_mes = "Unknown logger type [" + tokens[1] + "]";
    return false;
}

//----------------------------------------------------------------------
bool log_performance_by_command_str(
    const std::string & command,
    std::string & ret_mes)
{
    std::string separators = " ";
    std::vector< std::string > tokens;
    nv::index_common::Tokenizer::parse(command, separators, tokens);
    
    const mi::Size token_size = tokens.size();
    if((token_size <= 1) || (tokens[0] != "perf")){
        ret_mes = "Invalid command: perf [?|options]";
        return false;
    }

    if(tokens[1] == "?"){
        ret_mes = "Invalid command: perf [?|options]";
        INFO_LOG << "perf [options]\n"
                 << "    log the performance.\n"
                 << "    options\n"
                 << "    ?  help (this message)\n"
                 << "    frame=INT\n"
                 << "       measure the performance in INT frames.\n"
                 << "    system=bool\n"
                 << "       turn on|off the system logger. (default off)\n"
                 << "    per_host=bool\n"
                 << "       turn on|off the per_host logger. (default off)\n"
                 << "    per_span=bool\n"
                 << "       turn on|off the per_span logger. (default off)"
            ;
        return true;
    }

    mi::Sint32 nb_frame = 0;
    bool is_on[Performance_logger_base::LK_COUNT] = { false, false, false };

    for(mi::Size i = 1; i < token_size; ++i){
        const std::string key_val_sep = "=";
        std::vector< std::string > arg_tokens;
        nv::index_common::Tokenizer::parse(tokens[i], key_val_sep, arg_tokens);
        if(arg_tokens.size() != 2){
            ERROR_LOG << "perf: illegal args: " << tokens[i];
            continue;
        }

        if(arg_tokens[0] == "frame"){
            nb_frame = nv::index_common::get_sint32(arg_tokens[1]);
        }
        else if(arg_tokens[0] == "system"){
            is_on[Performance_logger_base::LK_System_global] = nv::index_common::get_bool(arg_tokens[1]);
        }
        else if(arg_tokens[0] == "per_host"){
            is_on[Performance_logger_base::LK_Per_host] = nv::index_common::get_bool(arg_tokens[1]);
        }
        else if(arg_tokens[0] == "per_span"){
            is_on[Performance_logger_base::LK_Per_span] = nv::index_common::get_bool(arg_tokens[1]);
        }
        else{
            ERROR_LOG << "perf: unrecognized key=value: " << tokens[i];
        }
    }

    bool is_success = true;
    if((1 < nb_frame) && (nb_frame < 10000)){
        // performance_log_state and initialization lock (peek_logger_set())
        mi::base::Lock::Block block(&Performance_logger_set_lock);

        INFO_LOG << "logging " << nb_frame << " frame performance.";
        Nvindex_AppData::instance()->set_remaining_performance_measuring_frames(nb_frame);
        
        Performance_logger_app_state & perf_log_state = Nvindex_AppData::instance()->peek_log_state();

        // turn on the loggers
        for(mi::Uint32 i = 0; i < Performance_logger_base::LK_COUNT; ++i){
            perf_log_state.set_logging_on(i, is_on[i]);
        }
        INFO_LOG << "performance log: system: " << is_on[Performance_logger_base::LK_System_global] 
                 << ", per_host: "              << is_on[Performance_logger_base::LK_Per_host] 
                 << ", per_span: "              << is_on[Performance_logger_base::LK_Per_span];
    }
    else{
        if(nb_frame == 1){
            // Special handling of the frame=1 case:
            //   1. set the logging on to GUI
            //   2a. Next frame, # of remaining frames = 0, turn of logging
            //   2b. GUI report the logging state changed as turn on logging usually next frame.
            // We don't distinguish the command is coming from GUI or command. This makes frame=1 
            // problem.
            ERROR_LOG << "perf: cannot handle one frame due to the GUI sync: please set frame more than 1.";
        }
        else{
            ERROR_LOG << "perf: illegal number of frames: " << nb_frame 
                      << ", the frame number is out of range. Ignored.";
        }
        is_success = false;
    }

    return is_success;
}

//----------------------------------------------------------------------
/// turn off all loggers
void turn_off_logging()
{
    // performance_log_state and initialization lock (peek_log_state())
    mi::base::Lock::Block block(&Performance_logger_set_lock);

    // turn off all the loggers
    Performance_logger_app_state & perf_log_state = Nvindex_AppData::instance()->peek_log_state();
    for(mi::Uint32 i = 0; i < Performance_logger_base::LK_COUNT; ++i){
        perf_log_state.set_logging_on(i, false);
    }
}

//----------------------------------------------------------------------
