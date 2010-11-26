
#ifndef __libpfetch_cmarshal_MARSHAL_H__
#define __libpfetch_cmarshal_MARSHAL_H__

#include	<glib-object.h>

G_BEGIN_DECLS

/* ENUM:STRING,POINTER,POINTER (libpfetch-cmarshal.list:1) */
extern void libpfetch_cmarshal_ENUM__STRING_POINTER_POINTER (GClosure     *closure,
                                                             GValue       *return_value,
                                                             guint         n_param_values,
                                                             const GValue *param_values,
                                                             gpointer      invocation_hint,
                                                             gpointer      marshal_data);

/* ENUM:VOID (libpfetch-cmarshal.list:2) */
extern void libpfetch_cmarshal_ENUM__VOID (GClosure     *closure,
                                           GValue       *return_value,
                                           guint         n_param_values,
                                           const GValue *param_values,
                                           gpointer      invocation_hint,
                                           gpointer      marshal_data);

G_END_DECLS

#endif /* __libpfetch_cmarshal_MARSHAL_H__ */

