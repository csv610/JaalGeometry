/****************************************************************************
** Meta object code from reading C++ file 'StructuredMeshDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/StructuredMeshDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'StructuredMeshDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JStructuredMeshDialog_t {
    QByteArrayData data[8];
    char stringdata0[82];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JStructuredMeshDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JStructuredMeshDialog_t qt_meta_stringdata_JStructuredMeshDialog = {
    {
QT_MOC_LITERAL(0, 0, 21), // "JStructuredMeshDialog"
QT_MOC_LITERAL(1, 22, 11), // "meshCreated"
QT_MOC_LITERAL(2, 34, 0), // ""
QT_MOC_LITERAL(3, 35, 7), // "genMesh"
QT_MOC_LITERAL(4, 43, 13), // "keyPressEvent"
QT_MOC_LITERAL(5, 57, 10), // "QKeyEvent*"
QT_MOC_LITERAL(6, 68, 1), // "e"
QT_MOC_LITERAL(7, 70, 11) // "closeDialog"

    },
    "JStructuredMeshDialog\0meshCreated\0\0"
    "genMesh\0keyPressEvent\0QKeyEvent*\0e\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JStructuredMeshDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   34,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       3,    0,   35,    2, 0x08 /* Private */,
       4,    1,   36,    2, 0x08 /* Private */,
       7,    0,   39,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 5,    6,
    QMetaType::Void,

       0        // eod
};

void JStructuredMeshDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JStructuredMeshDialog *_t = static_cast<JStructuredMeshDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->meshCreated(); break;
        case 1: _t->genMesh(); break;
        case 2: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        case 3: _t->closeDialog(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (JStructuredMeshDialog::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&JStructuredMeshDialog::meshCreated)) {
                *result = 0;
                return;
            }
        }
    }
}

const QMetaObject JStructuredMeshDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JStructuredMeshDialog.data,
      qt_meta_data_JStructuredMeshDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JStructuredMeshDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JStructuredMeshDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JStructuredMeshDialog.stringdata0))
        return static_cast<void*>(const_cast< JStructuredMeshDialog*>(this));
    if (!strcmp(_clname, "Ui::StructuredMeshDialog"))
        return static_cast< Ui::StructuredMeshDialog*>(const_cast< JStructuredMeshDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JStructuredMeshDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void JStructuredMeshDialog::meshCreated()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
